library(tidyverse)
library(stringr)

args = commandArgs(TRUE)
sample_name= args[1]
file_path= args[2]

#for designing primers flanking to the breakpoints
p3_out_file = paste0(file_path,"/gripss/",sample_name,"_p3_out")


p3_df = read_delim(p3_out_file, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
        mutate(is_blank = is.na(field)) %>%
        mutate(record_ordinal = cumsum(is_blank)) %>% 
        group_by(record_ordinal) %>%
        spread(field, value) %>%
        select(SEQUENCE_ID, PRIMER_LEFT_0_SEQUENCE, PRIMER_RIGHT_0_SEQUENCE) %>%
        mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
        mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
        mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
        filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)

#################################
#for designing primers overlapping the breakpoints
p3_out_file_overlap = paste0(file_path,"/gripss/",sample_name,"_p3_out_overlap")


p3_df_overlap = read_delim(p3_out_file_overlap, delim = "=", col_names = c("field", "value"), col_types = "cc") %>%
        mutate(is_blank = is.na(field)) %>%
        mutate(record_ordinal = cumsum(is_blank)) %>% 
        group_by(record_ordinal) %>%
        spread(field, value) %>%
        select(SEQUENCE_ID, PRIMER_LEFT_0_SEQUENCE, PRIMER_RIGHT_0_SEQUENCE) %>%
        mutate(PRIMER_LEFT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_LEFT","\n", PRIMER_LEFT_0_SEQUENCE,"\n")) %>%
        mutate(PRIMER_RIGHT_fasta = paste0(">", SEQUENCE_ID, " PRIMER_RIGHT","\n", PRIMER_RIGHT_0_SEQUENCE)) %>%
        mutate(PRIMER_PAIRS_fasta = paste0(PRIMER_LEFT_fasta, PRIMER_RIGHT_fasta)) %>%
        filter(is.na(PRIMER_LEFT_0_SEQUENCE) == FALSE)

###################
#bind both flanking and overlapping results
p3_bind = bind_rows(p3_df_overlap, p3_df) %>% ungroup() %>% 
        select(SEQUENCE_ID,PRIMER_PAIRS_fasta)


p3_bind = data.frame(separate(p3_bind, SEQUENCE_ID, c("QUAL","vcfID"), sep = "_"))

p3_bind = p3_bind[order(as.numeric(p3_bind$QUAL), decreasing = TRUE),]

writeLines(p3_bind$PRIMER_PAIRS_fasta, paste0(file_path,"/gripss/",sample_name,"_p3_out_full.fasta"))
#If both approach have primers in the same breakpoint, prioritising results from overlapping approach
p3_bind = p3_bind[!duplicated(p3_bind$vcfID),]

writeLines(head(p3_bind$PRIMER_PAIRS_fasta, 5), paste0(file_path,"/gripss/",sample_name,"_p3_out_top5.fasta"))
