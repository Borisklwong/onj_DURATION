library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(data.table)

for (sample_name in c("1-089_VD","3-080_ES","3-082_GT","3-084_JT","3-088_TN")) {




Cigar_metrics_folder="/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/FFPE_FF_BL/metrics/Cigar_metrics"
setwd(Cigar_metrics_folder)
file_name = "_merged_mark.bam.cigar_metrics"


#Plot bar charts to compare CIGAR=M between FF v.s FFPE (S) v.s buffy coat (BL)
for (cigar in c("M")) {

FF = as.data.frame(read.table(paste0(sample_name,"_FF",file_name), heade=TRUE, sep = "\t"))
BL = as.data.frame(read.table(paste0(sample_name,"_BL",file_name), header=TRUE, sep = "\t"))
S= as.data.frame(read.table(paste0(sample_name,"_S",file_name), header=TRUE, sep = "\t"))

FF_file = FF %>% filter(OPERATOR == cigar, LENGTH < 100) %>% mutate(SAMPLE = "FF")
BL_file = BL %>% filter(OPERATOR == cigar, LENGTH < 100) %>% mutate(SAMPLE = "BL")
S_file = S %>% filter(OPERATOR == cigar, LENGTH < 100) %>% mutate(SAMPLE = "S")

FULL_file = bind_rows(FF_file , BL_file, S_file) %>% mutate(ID=sample_name)

write.csv(FULL_file, paste0(sample_name,"_",cigar,"_cigar_metrics_FF_vs_FFPE.csv"))

plot = ggplot(FULL_file, aes(fill=SAMPLE, x=LENGTH, y=COUNT)) + 
        geom_col() + 
        ggtitle(paste0(sample_name," CIGAR=",cigar)) +
        scale_y_continuous(labels = comma)
pdf(paste0(sample_name,"_",cigar,"_plot.pdf"))
print(plot)
dev.off()

}

#Plot bar charts to compare CIGAR between FF v.s FFPE (S) v.s buffy coat (BL)
for (cigar in c("D","I","S","H")) {
  
  FF = as.data.frame(read.table(paste0(sample_name,"_FF",file_name), header=TRUE, sep = "\t"))
  BL = as.data.frame(read.table(paste0(sample_name,"_BL",file_name), header=TRUE, sep = "\t"))
  S= as.data.frame(read.table(paste0(sample_name,"_S",file_name), header=TRUE, sep = "\t"))
  
  FF_file = FF %>% filter(OPERATOR == cigar) %>% mutate(SAMPLE = "FF")
  BL_file = BL %>% filter(OPERATOR == cigar) %>% mutate(SAMPLE = "BL")
  S_file = S %>% filter(OPERATOR == cigar) %>% mutate(SAMPLE = "S")
  
  FULL_file = bind_rows(FF_file , BL_file, S_file)
  
  write.csv(FULL_file, paste0(sample_name,"_",cigar,"_cigar_metrics_FF_vs_FFPE.csv"))
  
  plot = ggplot(FULL_file, aes(fill=SAMPLE, x=LENGTH, y=COUNT)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle(paste0(sample_name," CIGAR=",cigar)) +
    scale_y_continuous(trans='log10', labels = comma)
  pdf(paste0(sample_name,"_",cigar,"_plot.pdf"))
  print(plot)
  dev.off()
  
}

#Plot line charts to compare insert_size_metrics between FF v.s FFPE (S) v.s buffy coat (BL)

insert_size_metrics="/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/FFPE_FF_BL/metrics/insert_size_metrics"
setwd(insert_size_metrics)
file_name = "_merged_mark.bam.insert_size_metrics"


  
FF = as.data.frame(read.table(paste0(sample_name,"_FF",file_name), header=TRUE, skip=10, sep = "\t"))
BL = as.data.frame(read.table(paste0(sample_name,"_BL",file_name), header=TRUE, skip=10, sep = "\t"))
S= as.data.frame(read.table(paste0(sample_name,"_S",file_name), header=TRUE, skip=10, sep = "\t"))

FF_file = FF %>% mutate(SAMPLE = "FF")
BL_file = BL %>% mutate(SAMPLE = "BL")
S_file = S %>% mutate(SAMPLE = "S")  

FULL_file = bind_rows(FF_file , BL_file, S_file)

write.csv(FULL_file, paste0(sample_name,"_insert_size_metrics_FF_vs_FFPE.csv"))

plot = ggplot(FULL_file, aes(group=SAMPLE, fill=SAMPLE, x=insert_size, y=All_Reads.fr_count)) + 
  geom_point(aes(color=SAMPLE)) + 
  geom_line(aes(linetype=SAMPLE, color =SAMPLE)) +
  ggtitle(paste0(sample_name,"_insert_size_metrics")) +
  scale_y_continuous(labels = comma)
  
##TODO: add peak value on the graph

pdf(paste0(sample_name,"_insert_size_plot.pdf"))
print(plot)
dev.off()

#Plot bar charts to compare tag_metrics between FF v.s FFPE (S) v.s buffy coat (BL)

tag_metrics="/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/FFPE_FF_BL/metrics/tag_metrics"
setwd(tag_metrics)
file_name = "_merged_mark.bam.tag_metrics"

FF = as.data.frame(read.table(paste0(sample_name,"_FF",file_name), header=TRUE, sep = "\t"))
BL = as.data.frame(read.table(paste0(sample_name,"_BL",file_name), header=TRUE, sep = "\t"))
S= as.data.frame(read.table(paste0(sample_name,"_S",file_name), header=TRUE,sep = "\t"))

FF_file = FF %>% mutate(SAMPLE = "FF", ID = sample_name)
BL_file = BL %>% mutate(SAMPLE = "BL", ID = sample_name)
S_file = S %>% mutate(SAMPLE = "S", ID = sample_name)  

FULL_file = bind_rows(FF_file , BL_file, S_file)

plot = ggplot(FULL_file, aes(fill=SAMPLE, x=TAG, y=COUNT, width=.75)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle(paste0(sample_name," tag_metrics")) +
  scale_y_continuous(trans='log10', labels = comma)

pdf(paste0(sample_name,"_tag_metrics_plot.pdf"))
print(plot)
dev.off()

#normalise count by NM
normalise_factor= FULL_file %>% filter(TAG == "NM")

normalise_factor = min(normalise_factor$COUNT)

FF_factor = FF %>% filter(TAG == "NM") %>% mutate(normaised_COUNT = normalise_factor/COUNT)
BL_factor = BL %>% filter(TAG == "NM") %>% mutate(normaised_COUNT = normalise_factor/COUNT)
S_factor = S %>% filter(TAG == "NM") %>% mutate(normaised_COUNT = normalise_factor/COUNT)

FF_file = FF_file %>% mutate(normaised_COUNT = COUNT*min(FF_factor$normaised_COUNT))
BL_file = BL_file %>% mutate(normaised_COUNT = COUNT*min(BL_factor$normaised_COUNT))
S_file = S_file %>% mutate(normaised_COUNT = COUNT*min(S_factor$normaised_COUNT))

FULL_file = bind_rows(FF_file , BL_file, S_file)

plot = ggplot(FULL_file, aes(fill=SAMPLE, x=TAG, y=normaised_COUNT, width=.75)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggtitle(paste0(sample_name," tag_metrics (normalised count)")) +
  scale_y_continuous(trans='log10', labels = comma)

pdf(paste0(sample_name,"_tag_metrics_plot_normalised.pdf"))
print(plot)
dev.off()

#Create table forFresh frozen vs FFPE count (FF vs FFPE)


FF_vs_FFPE_summarry = data.frame(
                      ID = sample_name,
                      TAG = FF_file$TAG, 
                      FF_COUNT=FF_file$COUNT, 
                      FFPE_COUNT=S_file$COUNT, 
                      "FF_vs_FFPE_(%)" =100*(S_file$COUNT/FF_file$COUNT),
                      FF_normaised=FF_file$normaised_COUNT,
                      FFPE_normaised=S_file$normaised_COUNT,
                      "FF_vs_FFPE_normaised(%)" =100*(S_file$normaised_COUNT/FF_file$normaised_COUNT)
                      )
write.csv(FF_vs_FFPE_summarry, paste0(sample_name,"_tag_metrics_FF_vs_FFPE.csv"))

}



merged_table = list.files(path=".", pattern="_tag_metrics_FF_vs_FFPE.csv$") %>% lapply(read_csv) %>% bind_rows
write.csv(merged_table, "tag_metrics_FF_vs_FFPE_summary.csv")

#Comparing tag across 5 samples
merged_table_tag = merged_table %>%
  melt(id.vars = c("ID", "TAG"),
       measure.vars = c("FF_normaised","FFPE_normaised"),
       variable.name = "SAMPLE",
       value.name = "normaised_COUNT")

plot = ggplot(merged_table_tag, aes(fill = SAMPLE, x=TAG, y=normaised_COUNT)) +
  geom_bar(position="dodge", stat="identity") +
  labs(
    title = "FF vs FFPE summary",
    subtitle = "(normalised COUNT)"
  ) +
  scale_y_continuous(trans='log10', labels = comma) +
  coord_flip() +
  facet_grid(rows = vars(ID), scales = "free_y", switch = "y", space = "free_y") +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    plot.title = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 7),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
  )
plot


pdf("FF_vs_FFPE_summary_plot.pdf")
print(plot)
dev.off()



