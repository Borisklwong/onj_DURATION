library(tidyverse)
library(stringr)


sample_name= "1-009_MC_S"
file_path= "/stornext/Home/data/allstaff/w/wong.b/onj_DURATION/processing/BGI_DURATION_processing/before_trimming/1-009_MC_S_process/1-009_MC_S_g-p-l_out"

#for designing primers flanking to the breakpoints
p3_out_file = paste0(file_path,"/gripss/",sample_name,"_p3_out")

load(paste0(sample_name,"_vcf_to_p3_df_new"))


p3_df = read_delim(p3_out_file, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
  mutate(is_blank = is.na(field)) %>%
  mutate(record_ordinal = cumsum(is_blank)) %>% 
  group_by(record_ordinal) %>%
  spread(field, value) %>%
  select(SEQUENCE_ID, SEQUENCE_TEMPLATE, PRIMER_LEFT_0_SEQUENCE, PRIMER_LEFT_0, PRIMER_LEFT_0_TM, PRIMER_RIGHT_0_SEQUENCE, PRIMER_RIGHT_0, PRIMER_RIGHT_0_TM, PRIMER_PAIR_0_LIBRARY_MISPRIMING, PRIMER_PAIR_0_PRODUCT_SIZE) %>%
  mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
  mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
  mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
  filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)

p3_df = data.frame(separate(p3_df, SEQUENCE_ID, c("QUAL","vcfId"), sep = "@"))
str_sub(p3_df$vcfId, -1, -1) = ""
p3_df = inner_join(p3_df, vcf_to_p3_df_new, by = "vcfId" , suffix = c("", ".df"))

#################################
#for designing primers overlapping the breakpoints
p3_out_file_overlap = paste0(file_path,"/gripss/",sample_name,"_p3_out_overlap")

load(paste0(sample_name,"_vcf_to_p3_df_overlap_new"))

p3_df_overlap = read_delim(p3_out_file_overlap, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
  mutate(is_blank = is.na(field)) %>%
  mutate(record_ordinal = cumsum(is_blank)) %>% 
  group_by(record_ordinal) %>%
  spread(field, value) %>%
  select(SEQUENCE_ID, SEQUENCE_TEMPLATE, PRIMER_LEFT_0_SEQUENCE, PRIMER_LEFT_0, PRIMER_LEFT_0_TM, PRIMER_RIGHT_0_SEQUENCE, PRIMER_RIGHT_0, PRIMER_RIGHT_0_TM, PRIMER_PAIR_0_LIBRARY_MISPRIMING, PRIMER_PAIR_0_PRODUCT_SIZE) %>%
  mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
  mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
  mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
  filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)
  
p3_df_overlap = data.frame(separate(p3_df_overlap, SEQUENCE_ID, c("QUAL","vcfId"), sep = "@"))
str_sub(p3_df_overlap$vcfId, -1, -1) = ""
p3_df_overlap = inner_join(p3_df_overlap, vcf_to_p3_df_overlap_new, by = "vcfId" , suffix = c("", ".df"))


###################
#summary for primers



#bind both flanking and overlapping results and extract amplicon sequence
p3_bind = bind_rows(p3_df, p3_df_overlap) %>% ungroup() %>%
  mutate(AMP_STR_POS = substr(p3_bind$PRIMER_LEFT_0,1,regexpr(",",p3_bind$PRIMER_LEFT_0)-1),
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


p3_bind = p3_bind[order(as.numeric(p3_bind$QUAL), decreasing = TRUE),] 


writeLines(p3_bind$PRIMER_PAIRS_fasta, paste0(file_path,"/gripss/",sample_name,"_p3_out_full.fasta"))
write.csv(p3_bind, paste0(file_path,"/gripss/",sample_name,"_p3_out_full.csv"))
writeLines(p3_bind$amplicon_fasta, paste0(file_path,"/gripss/",sample_name,"_ampicon_seq_full.fasta"))

#If both approach have primers in the same breakpoint, prioritising results from flanking approach
p3_bind = p3_bind[!duplicated(p3_bind$vcfId),]


writeLines(head(p3_bind$PRIMER_PAIRS_fasta, 5), paste0(file_path,"/gripss/",sample_name,"_p3_out_top5_new.fasta"))
write.csv(p3_bind, paste0(file_path,"/gripss/",sample_name,"_p3_out_top5_new.csv"))
writeLines(head(p3_bind$amplicon_fasta, 5), paste0(file_path,"/gripss/",sample_name,"_ampicon_seq_out_top5.fasta"))
