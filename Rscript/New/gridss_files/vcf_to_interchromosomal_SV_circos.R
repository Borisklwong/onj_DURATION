fileNames <- Sys.glob("*.vcf")

for (FF in fileNames) {

out_vcf.file=FF
out_bpgr = breakpointRanges(VariantAnnotation::readVcf(out_vcf.file))

#Only keep the interchromosomal SVs
out_bpgr_inter = out_bpgr[is.na(out_bpgr$svLen)]
#Pair the results
out_pairgr = breakpointgr2pairs(out_bpgr_inter)
# #sort by score
# out_pairgr_sort = out_pairgr[order(out_pairgr@first@elementMetadata@listData[["QUAL"]], decreasing = TRUE),]
# 
# out_pairgr_sort_name = paste(out_vcf.file,"_sort.csv")
# 
# 
# out_pairgr_sort_name_export = paste(out_vcf.file,".sv.bedpe")
# rtracklayer::export(out_pairgr_sort, con=out_pairgr_sort_name_export)


#Ploting circos plot
library(circlize)
library(ComplexHeatmap)
library(gridBase)
#change the format of the result from Granges into bedpe
out_bpgr_inter_bedpe = breakpointgr2bedpe(out_bpgr_inter)



plot_title=paste(gsub(pattern = ".gridss.somatic.vcf", "", basename(out_vcf.file)),".pdf",sep = "")
out_bpgr_inter_bedpe_filter = filter(out_bpgr_inter_bedpe, chrom1 %in% c(paste0("chr", 1:22), "chrX", "chrY", "M"))
out_bpgr_inter_bedpe_filter2 = filter(out_bpgr_inter_bedpe_filter, chrom2 %in% c(paste0("chr", 1:22), "chrX", "chrY", "M"))
filter_Q1000 = filter(out_bpgr_inter_bedpe_filter2, score >= 1000)
filter_Q500 = filter(out_bpgr_inter_bedpe_filter2, score <= 1000 & score >= 500)
filter_Qlow = filter(out_bpgr_inter_bedpe_filter2, score <= 500)
circle_size = unit(1, "snpc") # snpc unit gives you a square region
lgd_lines = Legend(at = c(">=1000", "500 to 999", ">499"), type = "lines", legend_gp = gpar(col = c("#D8441C", "#101110", "#28D81C"), lwd = 3), title_position = "topleft", title = "GRIDSS score", nrow = 1)


pdf(file=plot_title, onefile = FALSE)
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)

circos.initializeWithIdeogram()
title(gsub(pattern = ".gridss.somatic.vcf", "", basename(out_vcf.file)))
circos.genomicLink(filter_Q1000[,1:3], filter_Q1000[,4:6], col = add_transparency("#D8441C", transparency = 0.2))
circos.genomicLink(filter_Q500[,1:3], filter_Q500[,4:6], col = add_transparency("#101110", transparency = 0.4))
circos.genomicLink(filter_Qlow[,1:3], filter_Qlow[,4:6], col = add_transparency("#28D81C", transparency = 0.8))
upViewport()

draw(lgd_lines, y = unit(1, "npc") - circle_size, just = "bottom")


dev.off()

}