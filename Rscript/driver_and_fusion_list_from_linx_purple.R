library(dplyr)
library(readr)
library(plyr)
### WARNING: Do NOT run this script under LINX or PURPLE output folder #########

gpl_out_dir = ""
working_dir = ""
# Copy all the driver.catalog files to working dir
file.copy(Sys.glob(paste0(gpl_out_dir,"/*/*/*/*.driver.catalog.tsv")), working_dir)
# Copy all the .linx.fusion.tsv files to working dir
file.copy(Sys.glob(paste0(gpl_out_dir,"/*/*/*/*.linx.fusion.tsv")), working_dir)

# Summarising driver.catalog from purple
driver_files <- list.files(path = ".", pattern = ".driver.catalog.tsv")

for (i in 1:length(driver_files)) {
  
  
  sample = sub(pattern = ".driver.catalog.tsv", replacement = "\\1", basename(driver_files[i]))
  
  driver = as.data.frame(read.csv(driver_files[i], sep = "\t"))
  
  drive_list = data.frame(
    SampleID = rep(sample, nrow(driver)),
    chromosome = driver$chromosome,
    chromosomeBand = driver$chromosomeBand,
    Gene = driver$gene,
    DriverType = driver$driver,
    driverLikelihood = driver$driverLikelihood,
    category = driver$category
  )
  write.csv(drive_list, paste0(sample,"_driver_full.csv"))
}

driver_all = ldply(list.files(path = ".", pattern = "_driver_full.csv"), read.csv, header=TRUE)
write_csv(driver_all, "MESO_driver_gene_full.csv")


file.remove(list.files(path = ".", pattern = "_driver_full.csv"))

# Writing csv for the driver genes hits on the Gene List "gl"
driver_files <- list.files(path = ".", pattern = ".driver.catalog.tsv")

for (i in 1:length(driver_files)) {
  
  
  sample = sub(pattern = ".driver.catalog.tsv", replacement = "\\1", basename(driver_files[i]))
  
  driver = as.data.frame(read.csv(driver_files[i], sep = "\t"))
  # Please change the target gene list in gl
  gl = data.frame(gene = c("LATS1", "BAP1", "MST1", "NF2", "CDKN2A", "CDKN2B", "mTOR", "LASTS2", "STK3", "DDX3X", "DDX51", "SETD5", "SF3B1", "TRAF7", "SETDB1", "TP53", "TSC1", "TSC2", "ULK2", "SAV1", "TSC2", "ULK2", " SAV1"))
  
  gl_driver = semi_join(driver, gl, by = "gene")
  
  drive_list = data.frame(
    SampleID = rep(sample, nrow(gl_driver)),
    chromosome = gl_driver$chromosome,
    chromosomeBand = gl_driver$chromosomeBand,
    Gene = gl_driver$gene,
    DriverType = gl_driver$driver,
    driverLikelihood = gl_driver$driverLikelihood,
    category = gl_driver$category
  )
  write.csv(drive_list, paste0(sample,"_driver.csv"))
}

driver_all = ldply(list.files(path = ".", pattern = "_driver.csv"), read.csv, header=TRUE)
write_csv(driver_all, "MESO_driver_gene.csv")

file.remove(list.files(path = ".", pattern = "_driver.csv"))
file.remove(list.files(path = ".", pattern = ".driver.catalog.tsv"))
########################################
# Summarising linx.fusion from linx

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
write_csv(fus_all, "MESO_fusion_gene_full.csv")

file.remove(list.files(path = ".", pattern = "_fus_full.csv"))

# Writing csv for the fusion genes hits on the Gene List "gl"
fus_files <- list.files(path = ".", pattern = ".linx.fusion.tsv")

for (i in 1:length(fus_files)) {
  
  
  sample = sub(pattern = ".linx.fusion.tsv", replacement = "\\1", basename(fus_files[i]))
  
  fus = as.data.frame(read.csv(fus_files[i], sep = "\t"))
  # Please change the target gene list in gl
  gl = data.frame(GeneStart = c("LATS1", "BAP1", "MST1", "NF2", "CDKN2A", "CDKN2B", "mTOR", "LASTS2", "STK3", "DDX3X", "DDX51", "SETD5", "SF3B1", "TRAF7", "SETDB1", "TP53", "TSC1", "TSC2", "ULK2", "SAV1", "TSC2", "ULK2", " SAV1"),
                  GeneEnd = c("LATS1", "BAP1", "MST1", "NF2", "CDKN2A", "CDKN2B", "mTOR", "LASTS2", "STK3", "DDX3X", "DDX51", "SETD5", "SF3B1", "TRAF7", "SETDB1", "TP53", "TSC1", "TSC2", "ULK2", "SAV1", "TSC2", "ULK2", " SAV1"))
  
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
write_csv(fus_all, "MESO_fusion_gene.csv")

file.remove(list.files(path = ".", pattern = "_fus.csv"))
file.remove(list.files(path = ".", pattern = ".linx.fusion.tsv"))
