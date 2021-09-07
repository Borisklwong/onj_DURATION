#driver script: /home/users/allstaff/wong.b/onj_DURATION/scripts/gen_job_script/DURATION_circos_pairgr_autogen.sh

library(Rsamtools)
library(GenomicAlignments)
library(StructuralVariantAnnotation)
library(stringr)
library(dplyr)
library(data.table)

args = commandArgs(TRUE)
sample_name= args[1]
vcf_dir= args[2]
vcf_path= paste0(vcf_dir,"/gripss/",sample_name,"_49_PON_gridss_raw_with_hartwig_PON_new")
bam_path = paste0(vcf_dir,"/gridss/",sample_name,".assembly.bam.gridss.working")
# sample_name= "3-082_GT_S"
# vcf_path= "/stornext/Home/data/allstaff/w/wong.b/onj_DURATION/processing/BGI_DURATION_processing/3-082_GT_S_process/3-082_GT_S_g-p-l_out/gripss/3-082_GT_S_49_PON_gridss_raw_with_hartwig_PON_new"
# bam_path= "/stornext/Home/data/allstaff/w/wong.b/onj_DURATION/processing/BGI_DURATION_processing/3-082_GT_S_process/3-082_GT_S_g-p-l_out/gridss/3-082_GT_S.assembly.bam.gridss.working"

print(sample_name)
print(vcf_dir)
print(vcf_path)
print(bam_path)


out_vcf.file = paste0(vcf_path,"/",sample_name,"T.gripss.somatic.49.PON.filtered.vcf.gz")
#Filterout PON in vcf files and sort QUAL
vcf = VariantAnnotation::readVcf(out_vcf.file)
vcf_filt = vcf[filt(vcf) == "PASS"]
#vcf_filt_sort = vcf_filt[order(fixed(vcf_filt)$QUAL, decreasing = TRUE),]
#vcf_gr = rowRanges(vcf_filt_sort)

# Remove unpaired variants 
vcf_bpgr = breakpointRanges(vcf_filt)
#TO DO: filter variants  5kp
vcf_bpgr = vcf_bpgr[!(seqnames(vcf_bpgr) == seqnames(partner(vcf_bpgr)) & abs(start(vcf_bpgr) - start(partner(vcf_bpgr))) < 5000)]
#getting BEDID from vcf
vcf_bpgr$BEID = info(vcf_filt[vcf_bpgr$sourceId])$BEID
# #getting IHOMPOS from vcf
# vcf_bpgr$IHOMPOS = info(vcf_filt[vcf_bpgr$sourceId])$IHOMPOS
#getting RP from vcf
vcf_bpgr$RP = info(vcf_filt[vcf_bpgr$sourceId])$RP
#getting SR from vcf
vcf_bpgr$SR = info(vcf_filt[vcf_bpgr$sourceId])$SR



#Cross-match filtered vcf with GRIDSS assembly.bam
bam = paste0(bam_path,"/",sample_name,".assembly.bam.sv.bam")
bam = BamFile(bam)
params = ScanBamParam(which = resize(vcf_bpgr, 100, fix = "center"), what = scanBamWhat())
bam_filtered = scanBam(bam, param = params)


#as.data.frame(bam_filtered)

gal = readGAlignments(bam, use.names = TRUE, param = params)
gal_df = as.data.frame(gal)

##########################
# Create data.frame, from vcf files, matching bam and vcf
vcf_asm_df = data.frame(
  BEID = unlist(vcf_bpgr$BEID),
  vcfId = rep(vcf_bpgr$sourceId, elementNROWS(vcf_bpgr$BEID))
) %>% 
  inner_join(as.data.frame(vcf_bpgr), by = c("vcfId" = "sourceId"), suffix = c("", ".vcf")) %>%
  inner_join(gal_df, by = c("BEID" = "qname"), suffix = c(".vcf", ".asm")) %>%
  mutate(asm_break_pos = ifelse(strand.vcf == "+", end.asm, start.asm)) %>%
  filter(asm_break_pos >= start.vcf & asm_break_pos <= end.vcf & seqnames.vcf == seqnames.asm) %>%
  distinct() %>%
  # Remove calls with > 1 "S"s in cigar(i.e. more than 1 "S"s in cigar)
  filter(str_count(cigar, "S") == 1) %>%
  # Only keep insertion length <= 3 bps
  filter(insLen <= 3) %>%
  filter(RP >= 1) %>%
  filter(SR >= 2) %>%
  filter(HOMLEN <=3)

#only keep calls with svLen >= 10000 bps or interchromosomal
vcf_asm_df = vcf_asm_df[!vcf_asm_df$svLen %between% c(-10000,10000) | is.na(vcf_asm_df$svLen),]

#sort vcf_asm_df by GRIDSS score
vcf_asm_df = vcf_asm_df[order(vcf_asm_df$QUAL, decreasing = TRUE),]

vcf_asm_df_partner = data.frame(
  partner_vcfId = vcf_asm_df$partner,
  partner_seqnames = vcf_asm_df$seqnames.vcf
)  

vcf_asm_df = inner_join(vcf_asm_df, vcf_asm_df_partner, by = c("vcfId" = "partner_vcfId"))

# Identify break point position for primer 3
#length of the local or distal assembly after including "I" and excluding "D" in cigar (i.e. after soft clipping)
qlen = cigarWidthAlongQuerySpace(vcf_asm_df$cigar, after.soft.clipping = TRUE)

p3_seq_target_start = ifelse(vcf_asm_df$strand.vcf == "+", qlen, vcf_asm_df$qwidth - qlen - vcf_asm_df$insLen)


#Remove the seq 
# Make this as a function, use  test_that to check
p3_seq_target =  ifelse(vcf_asm_df$insLen == 0, paste0(p3_seq_target_start,", 1"), paste0(p3_seq_target_start,",",vcf_asm_df$insLen))

vcf_asm_df = vcf_asm_df %>% mutate(P3_SEQ_TARGET = p3_seq_target, P3_TARGET = p3_seq_target_start)



write.csv(vcf_asm_df, paste0(vcf_path,"/",sample_name,"_vcf_asm_df.csv"))

# Output result for Primer3 input format, primers flanking to breakpoints
vcf_to_p3_df = vcf_asm_df %>% mutate(p3_record = paste0(
  "SEQUENCE_ID=",vcf_asm_df$QUAL,"@",vcf_asm_df$vcfId,"@",sample_name,"@chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames,"\n",
  "SEQUENCE_TEMPLATE=",vcf_asm_df$seq,"\n",
  "SEQUENCE_TARGET=",vcf_asm_df$P3_SEQ_TARGET,"\n",
  "="),
  SEQ_NAME = paste0(vcf_asm_df$QUAL,"@",vcf_asm_df$vcfId,"@",sample_name,"@chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames))
#remove the "o" or "h" after vcfId
str_sub(vcf_to_p3_df$vcfId, -1, -1) = ""
vcf_to_p3_df_new = vcf_to_p3_df[!duplicated(vcf_to_p3_df$vcfId),]
save(vcf_to_p3_df_new, file = paste0(sample_name,"_vcf_to_p3_df_new"))
writeLines(vcf_to_p3_df_new$p3_record, paste0(vcf_path,"/",sample_name,"_vcf_to_p3.txt"))

# Output result for Primer3 input format, primers overlapping to breakpoints
vcf_to_p3_df_overlap = vcf_asm_df %>% mutate(p3_record = paste0(
  "SEQUENCE_ID=",vcf_asm_df$QUAL,"@",vcf_asm_df$vcfId,"@",sample_name,"@chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames,"@overlap","\n",
  "SEQUENCE_TEMPLATE=",vcf_asm_df$seq,"\n",
  "SEQUENCE_OVERLAP_JUNCTION_LIST=",vcf_asm_df$P3_TARGET,"\n",
  "="),
  SEQ_NAME = paste0(vcf_asm_df$QUAL,"@",vcf_asm_df$vcfId,"@",sample_name,"@chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames,"@overlap"))
#remove the "o" or "h" after vcfId
str_sub(vcf_to_p3_df_overlap$vcfId, -1, -1) = ""
vcf_to_p3_df_overlap_new = vcf_to_p3_df_overlap[!duplicated(vcf_to_p3_df_overlap$vcfId),]
save(vcf_to_p3_df_overlap_new, file = paste0(sample_name,"_vcf_to_p3_df_overlap_new"))
writeLines(vcf_to_p3_df_overlap_new$p3_record, paste0(vcf_path,"/",sample_name,"_vcf_to_p3_overlap.txt"))