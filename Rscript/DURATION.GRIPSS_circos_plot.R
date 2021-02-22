#driver script: /home/users/allstaff/wong.b/onj_DURATION/scripts/gen_job_script/DURATION_circos_pairgr_autogen.sh

args = commandArgs(TRUE)
sample_name= args[1]
file_path= args[2]
vcf_path = paste0(file_path,"/gripss")
library(StructuralVariantAnnotation)
library(dplyr)

out_vcf.file = paste0(vcf_path,"/",sample_name,"T.gripss.somatic.filtered.vcf.gz")
out_bpgr = breakpointRanges(VariantAnnotation::readVcf(out_vcf.file))
#Filter the results in PON
out_bpgr_PON = out_bpgr[out_bpgr$FILTER %in% "PASS"]
#Only keep the interchromosomal SVs
out_bpgr_inter = out_bpgr_PON[is.na(out_bpgr_PON$svLen)]
#Pair the results
out_pairgr = breakpointgr2pairs(out_bpgr_inter)
#sort by score
out_pairgr_sort = out_pairgr[order(out_pairgr@first@elementMetadata@listData[["QUAL"]], decreasing = TRUE),]
write.csv(out_pairgr_sort, paste0(vcf_path,"/",sample_name,"_out_pairgr_sort.csv"))

#Ploting circos plot
library(circlize)
library(ComplexHeatmap)
library(gridBase)
#change the format of the result from Granges into bedpe
out_bpgr_inter_bedpe = breakpointgr2bedpe(out_bpgr_inter)
out_bpgr_inter_bedpe_sort = out_bpgr_inter_bedpe[order(out_bpgr_inter_bedpe$score, decreasing = TRUE),]
write.csv(out_bpgr_inter_bedpe_sort, paste0(vcf_path,"/",sample_name,"_out_bpgr_inter_bedpe.csv"))
  
plot_title=paste(gsub(pattern = ".repeatmasker.somatic.filtered.vcf", "", sample_name),".pdf",sep = "")
out_bpgr_inter_bedpe_filter = filter(out_bpgr_inter_bedpe, chrom1 %in% c(paste0(1:22), "X", "Y", "M"))
out_bpgr_inter_bedpe_filter2 = filter(out_bpgr_inter_bedpe_filter, chrom2 %in% c(paste0(1:22), "X", "Y", "M"))
#Changing Chromosome Notation, Adding "chr" in front of chrom 1 and chrom2
out_bpgr_inter_bedpe_filter2$chrom1 = paste("chr",out_bpgr_inter_bedpe_filter2$chrom1, sep="")
out_bpgr_inter_bedpe_filter2$chrom2 = paste("chr",out_bpgr_inter_bedpe_filter2$chrom2, sep="")

filter_Q1000 = filter(out_bpgr_inter_bedpe_filter2, score >= 1000)
filter_Q500 = filter(out_bpgr_inter_bedpe_filter2, score <= 1000 & score >= 500)
filter_Qlow = filter(out_bpgr_inter_bedpe_filter2, score <= 500)
circle_size = unit(1, "snpc") # snpc unit gives you a square region
lgd_lines = Legend(at = c(">=1000", "500 to 999", ">499"), type = "lines", legend_gp = gpar(col = c("#D8441C", "#101110", "#28D81C"), lwd = 3), title_position = "topleft", title = "GRIDSS score", nrow = 1)
  
  
pdf(file=plot_title, onefile = FALSE)
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)
  
circos.initializeWithIdeogram()
title(gsub(pattern = ".repeatmasker.somatic.filtered.vcf", "", sample_name))
circos.genomicLink(filter_Q1000[,1:3], filter_Q1000[,4:6], col = add_transparency("#D8441C", transparency = 0.1))
circos.genomicLink(filter_Q500[,1:3], filter_Q500[,4:6], col = add_transparency("#101110", transparency = 0.9))
circos.genomicLink(filter_Qlow[,1:3], filter_Qlow[,4:6], col = add_transparency("#28D81C", transparency = 0.8))
upViewport()
  
draw(lgd_lines, y = unit(1, "npc") - circle_size, just = "bottom")
  
  
dev.off()
  

