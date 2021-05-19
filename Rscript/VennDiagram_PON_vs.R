library("VennDiagram")

args = commandArgs(TRUE)
sample_name= args[1]
file_path= args[2]

hartwig_PON = paste0(file_path,"/",sample_name,"_g-p-l_out","/gripss","/",sample_name,"_out_bpgr_inter_bedpe.csv")
PON_49 = paste0(file_path,"/",sample_name,"_g-p-l_out","/gripss","/",sample_name,"_49_PON","/",sample_name,"_out_bpgr_inter_bedpe.csv")

hartwig_PON=subset(read.csv(hartwig_PON), select = -c(X))
PON_49=subset(read.csv(PON_49), select = -c(X))

venn.diagram(
  x = list(hartwig_PON$name, PON_49$name),
  category.names = c("hartwig_PON","PON_49_BGI"),
  filename = paste0(sample_name,"_Venn.png"),
  output = TRUE,
  fill = c("pink","yellow")
)