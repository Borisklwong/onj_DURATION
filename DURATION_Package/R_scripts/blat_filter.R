library("Biostrings")
library(dplyr)

args = commandArgs(TRUE)
sample_name= args[1]
vcf_dir= args[2]
file_path=paste0(vcf_dir,"/gripss/",sample_name,"_49_PON_gridss_raw_with_hartwig_PON_new")
# sample_name="3-082_GT_S"
# 
# file_path="."

psl = paste0(file_path,"/",sample_name,"_p3_out_full.psl")
csv = as.data.frame(read.csv(paste0(file_path,"/",sample_name,"_p3_out_full.csv")))

blat_df = as.data.frame(read.table(psl, header = FALSE, sep = "\t", skip =5)) %>%
  rename("SEQ_NAME" = "V10")

blat_count = as.data.frame(table(blat_df$SEQ_NAME)) %>%
  rename("PIRMER_NAME" = "Var1")

blat_count$SEQ_NAME = gsub("_PRIMER_LEFT","", as.character(blat_count$PIRMER_NAME))
blat_count$SEQ_NAME = gsub("_PRIMER_RIGHT","", as.character(blat_count$SEQ_NAME))

final_csv = full_join(csv, blat_count, by = "SEQ_NAME")

fastaFile = readDNAStringSet(paste0(file_path,"/",sample_name,"_p3_out_full.fasta"))
seq_name = names(fastaFile)
sequence = paste(fastaFile)
fasta_df = data.frame(seq_name, sequence) %>%
  rename("PIRMER_NAME" ="seq_name")


df_join = inner_join(fasta_df, blat_count, by = "PIRMER_NAME")


final_df = inner_join(final_csv, df_join, by = "PIRMER_NAME")

write.csv(final_df, paste0(file_path,"/",sample_name,"_p3_out_full_counted.csv"))
