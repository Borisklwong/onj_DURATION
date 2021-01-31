library(Rsamtools)
library(GenomicAlignments)
library(VariantAnnotation)
library(splitstackshape)
library(splitstackshape)

#Filterout PON in vcf files and sort QUAL
vcf_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/1-009_MC_S_process/1-009_MC_S_g-p-l_out/gripss"
vcf = paste0(vcf_dir,"/","1-009_MC_ST.gripss.somatic.filtered.vcf.gz")
vcf = VariantAnnotation::readVcf(vcf)
vcf_filt = vcf[filt(vcf) == "PASS"]
vcf_filt_sort = vcf_filt[order(fixed(vcf_filt)$QUAL, decreasing = TRUE),]
#vcf_gr = rowRanges(vcf_filt_sort)

# Remove unpaired variants 
vcf_bpgr = breakpointRanges(vcf_filt_sort)

#Cross-match filtered vcf with GRIDSS assembly.bam
bam_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/1-009_MC_S_process/1-009_MC_S_g-p-l_out/gridss/1-009_MC_S.assembly.bam.gridss.working"
bam = paste0(bam_dir,"/","1-009_MC_S.assembly.bam.sv.bam")
bam = BamFile(bam)
params = ScanBamParam(which = resize(vcf_bpgr, 100, fix = "center"), what = scanBamWhat())
bam_filtered = scanBam(bam, param = params)


#as.data.frame(bam_filtered)

gal = readGAlignments(bam, use.names = TRUE, param = params)
gal_df = as.data.frame(gal)

##########################
# Create data.frame, from vcf files, before remove unpaired variants
vcf_df = data.frame(
  vcfId = names(vcf_filt_sort),
  BEID = unstrsplit(info(vcf_filt_sort)$BEID, sep = ","),
  QUAL = fixed(vcf_filt_sort)$QUAL,
  
  ALT = fixed(vcf_filt_sort)$ALT
  )
# Create data.frame, from breakpointRanges
vcf_paried_df = data.frame(
  vcfId = vcf_bpgr$sourceId,
  svLen = vcf_bpgr$svLen,
  insSeq = vcf_bpgr$insSeq,
  insLen = vcf_bpgr$insLen,
  HOMLEN = vcf_bpgr$HOMLEN
  )
# Combine Date from two data frame, Retain only paired variants
combine_df = dplyr::inner_join(vcf_df, vcf_paried_df, by = "vcfId")
slipt_df = cSplit(combine_df, "BEID", ",")


slipt_1 = rename(slipt_df, qname = BEID_1)
join_df_1 = inner_join(slipt_1, gal_df, by = "qname")
join_df_1 = rename(join_df_1, BEID_1 = qname)

slipt_2 = rename(slipt_df, qname = BEID_2)
join_df_2 = inner_join(slipt_2, gal_df, by = "qname")
join_df_2 = rename(join_df_2, BEID_2 = qname)

slipt_3 = rename(slipt_df, qname = BEID_3)
join_df_3 = inner_join(slipt_3, gal_df, by = "qname")
join_df_3 = rename(join_df_3, BEID_3 = qname)

join_df_1_2_3 = join_df_1 %>% union(join_df_2) %>% union(join_df_3) %>% arrange(desc(QUAL))

final_df = data.frame(
  vcfId = join_df_1_2_3$vcfId,
  QUAL = join_df_1_2_3$QUAL,
  ALT = join_df_1_2_3$ALT.value,
  insSeq= join_df_1_2_3$insSeq,
  insLen = join_df_1_2_3$insLen,
  svLen = join_df_1_2_3$svLen,
  HOMLEN = join_df_1_2_3$HOMLEN,
  strand = join_df_1_2_3$strand,
  cigar= join_df_1_2_3$cigar.1,
  seq = join_df_1_2_3$seq
)
## TO DO: remove duplicate rows