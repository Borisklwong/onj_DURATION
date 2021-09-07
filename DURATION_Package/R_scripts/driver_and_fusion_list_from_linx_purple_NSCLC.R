library(tidyverse)
genes_of_interest = c("EGFR", "KRAS", "ALK", "ROS1", "BRAF", "NTRK", "MET", "RET")
gpl_out_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/FFPE_only"
working_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/onj_DURATION/processing/BGI_DURATION_processing/FFPE_only/FFPE_only_drivers"
project_name = "DURATION_FFPE"


# Summarising driver.catalog from purple

read_driver_catalog = function (filename) {
  sample = sub(pattern = ".driver.catalog.tsv", replacement = "\\1", basename(filename))
  drive_list = read_tsv(filename) %>% 
    dplyr::select(
      chromosome, 
      chromosomeBand, 
      gene, 
      driver, 
      driverLikelihood, 
      category) %>% 
    mutate(SampleID=sample)
  return(drive_list)
}

driver_all = bind_rows(lapply(list.files(path = gpl_out_dir, pattern = ".driver.catalog.tsv", recursive = TRUE, full.names = TRUE), read_driver_catalog))
write_csv(driver_all, paste0(project_name,"_driver_gene_full.csv"))

driver_filtered = driver_all %>% 
  filter(gene %in% genes_of_interest)

########################################
# Summarising linx.fusion from linx
#######TO DO: change the code below ########

fus_files <- list.files(path = ".", pattern = ".linx.fusion.tsv")

for (i in 1:length(fus_files)) {
  
  
  sample = sub(pattern = ".linx.fusion.tsv", replacement = "\\1", basename(fus_files[i]))
  
  fus = as.data.frame(read.csv(fus_files[i], sep = "\t"))
  
  fus_list = data.frame(
    SampleID = rep(sample, nrow(fus)),
    FusionName = fus$Name,
    GeneStart = fus$GeneStart,
    StartExon = fus$GeneContextStart,
    GeneEnd = fus$GeneEnd,
    EndExon = fus$GeneContextEnd,
    JunctionCopyNumber = fus$JunctionCopyNumber
  )
  
  write.csv(fus_list, paste0(sample,"_fus_full.csv"))
}

fus_all = ldply(list.files(path = ".", pattern = "_fus_full.csv"), read.csv, header=TRUE)
write_csv(fus_all, paste0(project_name,"_fusion_gene_full.csv"))

file.remove(list.files(path = ".", pattern = "_fus_full.csv"))

# Writing csv for the fusion genes hits on the Gene List "gl"
fus_files <- list.files(path = ".", pattern = ".linx.fusion.tsv")

for (i in 1:length(fus_files)) {
  
  
  sample = sub(pattern = ".linx.fusion.tsv", replacement = "\\1", basename(fus_files[i]))
  
  fus = as.data.frame(read.csv(fus_files[i], sep = "\t"))
  # Please change the target gene list in gl
  gl = data.frame(GeneStart = genes_of_interest,
                  GeneEnd = genes_of_interest)
  
  gl_fus_start = semi_join(fus, gl, by = "GeneStart")
  gl_fus_end = semi_join(fus, gl, by = "GeneEnd")
  gl_fus = union(gl_fus_start, gl_fus_end)
  
  fus_list = data.frame(
    SampleID = rep(sample, nrow(gl_fus)),
    FusionName = gl_fus$Name,
    GeneStart = gl_fus$GeneStart,
    StartExon = gl_fus$GeneContextStart,
    GeneEnd = gl_fus$GeneEnd,
    EndExon = gl_fus$GeneContextEnd,
    JunctionCopyNumber = gl_fus$JunctionCopyNumber
  )
  
  write.csv(fus_list, paste0(sample,"_fus.csv"))
}

fus_all = ldply(list.files(path = ".", pattern = "_fus.csv"), read.csv, header=TRUE)
write_csv(fus_all, paste0(sample,"_fusion_gene.csv"))

file.remove(list.files(path = ".", pattern = "_fus.csv"))
file.remove(list.files(path = ".", pattern = ".linx.fusion.tsv"))
