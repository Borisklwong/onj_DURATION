library(tidyverse)
library(stringr)
library(stringi)

args = commandArgs(TRUE)
sample_name= args[1]
vcf_dir= args[2]
R_temp_file_dir = args[3]
file_path=paste0(vcf_dir,"/gripss/",sample_name,"_49_PON_gridss_raw_with_hartwig_PON_new")
# sample_name= "3-082_GT_S"
# file_path= "/stornext/Home/data/allstaff/w/wong.b/onj_DURATION/processing/BGI_DURATION_processing/3-082_GT_S_process/3-082_GT_S_g-p-l_out/gripss/3-082_GT_S_49_PON_gridss_raw_with_hartwig_PON_new"
# R_temp_file_dir = "/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/post_gridss_FFPE_only_primers"

print(sample_name)
print(R_temp_file_dir)
print(file_path)



#for designing primers flanking to the breakpoints
p3_out_file = paste0(file_path,"/",sample_name,"_p3_out")
print(p3_out_file)
load(paste0(R_temp_file_dir,"/",sample_name,"_vcf_to_p3_df_new"))

# #get the abs sum of IHOMPOS
# IHOMPOS_sum = as.data.frame(unlist(lapply(lapply(vcf_to_p3_df_new$IHOMPOS, abs), sum)))
# colnames(IHOMPOS_sum) = "IHOMPOS"
# vcf_to_p3_df_new$IHOMPOS = IHOMPOS_sum


p3_df = read_delim(p3_out_file, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
  mutate(is_blank = is.na(field)) %>%
  mutate(record_ordinal = cumsum(is_blank)) %>% 
  group_by(record_ordinal) %>%
  spread(field, value) %>%
  dplyr::select(SEQUENCE_ID, SEQUENCE_TEMPLATE, PRIMER_LEFT_0_SEQUENCE, PRIMER_LEFT_0, PRIMER_LEFT_0_TM, PRIMER_RIGHT_0_SEQUENCE, PRIMER_RIGHT_0, PRIMER_RIGHT_0_TM, PRIMER_PAIR_0_LIBRARY_MISPRIMING, PRIMER_PAIR_0_PRODUCT_SIZE) %>%
  mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, "_PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
  mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, "_PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
  mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
  filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)

p3_df = data.frame(separate(p3_df, SEQUENCE_ID, c("QUAL","vcfId"), sep = "@"))
str_sub(p3_df$vcfId, -1, -1) = ""
p3_df = inner_join(p3_df, vcf_to_p3_df_new, by = "vcfId" , suffix = c("", ".df"))
# p3_df = subset(p3_df, select = -c(IHOMPOS))



#################################
#for designing primers overlapping the breakpoints
p3_out_file_overlap = paste0(file_path,"/",sample_name,"_p3_out_overlap")

load(paste0(R_temp_file_dir,"/",sample_name,"_vcf_to_p3_df_overlap_new"))

# #get the abs sum of IHOMPOS
# IHOMPOS_sum_overlap = as.data.frame(unlist(lapply(lapply(vcf_to_p3_df_overlap_new$IHOMPOS, abs), sum)))
# colnames(IHOMPOS_sum_overlap) = "IHOMPOS"
# vcf_to_p3_df_overlap_new$IHOMPOS = IHOMPOS_sum

p3_df_overlap = read_delim(p3_out_file_overlap, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
  mutate(is_blank = is.na(field)) %>%
  mutate(record_ordinal = cumsum(is_blank)) %>% 
  group_by(record_ordinal) %>%
  spread(field, value) %>%
  dplyr::select(SEQUENCE_ID, SEQUENCE_TEMPLATE, PRIMER_LEFT_0_SEQUENCE, PRIMER_LEFT_0, PRIMER_LEFT_0_TM, PRIMER_RIGHT_0_SEQUENCE, PRIMER_RIGHT_0, PRIMER_RIGHT_0_TM, PRIMER_PAIR_0_LIBRARY_MISPRIMING, PRIMER_PAIR_0_PRODUCT_SIZE) %>%
  mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, "_PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
  mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, "_PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
  mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
  filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)

p3_df_overlap = data.frame(separate(p3_df_overlap, SEQUENCE_ID, c("QUAL","vcfId"), sep = "@"))
str_sub(p3_df_overlap$vcfId, -1, -1) = ""
p3_df_overlap = inner_join(p3_df_overlap, vcf_to_p3_df_overlap_new, by = "vcfId" , suffix = c("", ".df"))
# p3_df_overlap = subset(p3_df_overlap, select = -c(IHOMPOS))

###################
#summary for primers



#bind both flanking and overlapping results and extract amplicon sequence
p3_bind = bind_rows(p3_df, p3_df_overlap) %>% ungroup()

p3_bind = mutate(p3_bind, AMP_STR_POS = substr(p3_bind$PRIMER_LEFT_0,1,regexpr(",",p3_bind$PRIMER_LEFT_0)-1),
                 AMP_END_POS = substr(p3_bind$PRIMER_RIGHT_0,1,regexpr(",",p3_bind$PRIMER_RIGHT_0)-1))
p3_bind$AMP_STR_POS = as.numeric(p3_bind$AMP_STR_POS)+1
p3_bind$AMP_END_POS = as.numeric(p3_bind$AMP_END_POS)+1
p3_bind = mutate(p3_bind, amplicon_seq = substr(p3_bind$SEQUENCE_TEMPLATE, p3_bind$AMP_STR_POS, p3_bind$AMP_END_POS)) %>%
  mutate(amplicon_fasta = paste0(">",p3_bind$SEQ_NAME,"\n", amplicon_seq))
#######
start = p3_bind$AMP_STR_POS
start_len = nchar(p3_bind$PRIMER_LEFT_0_SEQUENCE)
end = p3_bind$AMP_END_POS
end_len = nchar(p3_bind$PRIMER_RIGHT_0_SEQUENCE)
break_pos = p3_bind$P3_TARGET


before_F_primer = substr(p3_bind$SEQUENCE_TEMPLATE,0,start-1)
F_primer = tolower(substr(p3_bind$SEQUENCE_TEMPLATE,start,start-1+start_len))
R_primer= tolower(substr(p3_bind$SEQUENCE_TEMPLATE,end+1-end_len,end))
berween_primer = substr(p3_bind$SEQUENCE_TEMPLATE,start+start_len,end-end_len)
after_R_primer = substr(p3_bind$SEQUENCE_TEMPLATE,end+1,nchar(p3_bind$SEQUENCE_TEMPLATE))


visaulise_seq = paste0(before_F_primer,F_primer,berween_primer,R_primer,after_R_primer)
stri_sub(visaulise_seq, break_pos+1, break_pos) <- "__"

p3_bind = mutate(p3_bind, visaulise_seq = visaulise_seq)
##################


p3_bind = p3_bind[order(as.numeric(p3_bind$QUAL), decreasing = TRUE),] %>% subset(select = -c(BEID.vcf,P3_SEQ_TARGET))


writeLines(p3_bind$PRIMER_PAIRS_fasta, paste0(file_path,"/",sample_name,"_p3_out_full.fasta"))
write.csv(p3_bind, paste0(file_path,"/",sample_name,"_p3_out_full.csv"))
writeLines(p3_bind$amplicon_fasta, paste0(file_path,"/",sample_name,"_ampicon_seq_full.fasta"))

#If both approach have primers in the same breakpoint, prioritising results from flanking approach
p3_bind = p3_bind[!duplicated(p3_bind$vcfId),]


writeLines(head(p3_bind$PRIMER_PAIRS_fasta, 10), paste0(file_path,"/",sample_name,"_p3_out_top10_new.fasta"))
write.csv(p3_bind, paste0(file_path,"/",sample_name,"_p3_out_top5_new.csv"))
writeLines(head(p3_bind$amplicon_fasta, 10), paste0(file_path,"/",sample_name,"_ampicon_seq_out_top10.fasta"))
