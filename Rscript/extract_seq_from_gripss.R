#driver script: /home/users/allstaff/wong.b/onj_DURATION/scripts/gen_job_script/DURATION_circos_pairgr_autogen.sh

library(Rsamtools)
library(GenomicAlignments)
library(StructuralVariantAnnotation)
library(stringr)
library(dplyr)

args = commandArgs(TRUE)
sample_name= args[1]
file_path= args[2]
vcf_path = paste0(file_path,"/gripss")
bam_path = paste0(file_path,"/gridss/",sample_name,".assembly.bam.gridss.working")

out_vcf.file = paste0(vcf_path,"/",sample_name,"T.gripss.somatic.filtered.vcf.gz")
#Filterout PON in vcf files and sort QUAL
vcf = VariantAnnotation::readVcf(out_vcf.file)
vcf_filt = vcf[filt(vcf) == "PASS"]
#vcf_filt_sort = vcf_filt[order(fixed(vcf_filt)$QUAL, decreasing = TRUE),]
#vcf_gr = rowRanges(vcf_filt_sort)

# Remove unpaired variants 
vcf_bpgr = breakpointRanges(vcf_filt)
#TO DO: filter variants  5kp
vcf_bpgr = vcf_bpgr[!(seqnames(vcf_bpgr) == seqnames(partner(vcf_bpgr)) & abs(start(vcf_bpgr) - start(partner(vcf_bpgr))) < 5000)]
vcf_bpgr$BEID = info(vcf_filt[vcf_bpgr$sourceId])$BEID



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
# Only keep interchromosomal calls
  filter(is.na(svLen) == TRUE)

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
  "SEQUENCE_ID=",vcf_asm_df$vcfId,"_",sample_name,"_",vcf_asm_df$QUAL,"_chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames,"\n",
  "SEQUENCE_TEMPLATE=",vcf_asm_df$seq,"\n",
  "SEQUENCE_TARGET=",vcf_asm_df$P3_SEQ_TARGET,"\n",
  "="))
str_sub(vcf_to_p3_df$vcfId, -1, -1) = ""
vcf_to_p3_df_new = vcf_to_p3_df[!duplicated(vcf_to_p3_df$vcfId),]

writeLines(vcf_to_p3_df_new$p3_record, paste0(vcf_path,"/",sample_name,"_vcf_to_p3.txt"))

# Output result for Primer3 input format, primers overlapping to breakpoints
vcf_to_p3_df_overlap = vcf_asm_df %>% mutate(p3_record = paste0(
  "SEQUENCE_ID=",vcf_asm_df$vcfId,"_",sample_name,"_",vcf_asm_df$QUAL,"_chr",vcf_asm_df$seqnames.vcf,"-",vcf_asm_df$partner_seqnames,"_overlap","\n",
  "SEQUENCE_TEMPLATE=",vcf_asm_df$seq,"\n",
  "SEQUENCE_OVERLAP_JUNCTION_LIST=",vcf_asm_df$P3_TARGET,"\n",
  "="))
str_sub(vcf_to_p3_df_overlap$vcfId, -1, -1) = ""
vcf_to_p3_df_overlap_new = vcf_to_p3_df_overlap[!duplicated(vcf_to_p3_df_overlap$vcfId),]

writeLines(vcf_to_p3_df_overlap_new$p3_record, paste0(vcf_path,"/",sample_name,"_vcf_to_p3_overlap.txt"))