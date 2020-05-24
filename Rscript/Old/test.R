# example sciprt performing simple filtering of variants to somatic calls
#source("https://bioconductor.org/biocLite.R")
BiocManager::install("VariantAnnotation")
library(stringr)
library(VariantAnnotation)
library(devtools)
install_github("PapenfussLab/StructuralVariantAnnotation")
library(StructuralVariantAnnotation)

FR07886867_vcf_example <- readVcf(FR07886867_vcf)
# filter out low quality calls
#vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
# somatic calls have no support in the normal
#somatic_vcf <- vcf[geno(vcf)$QUAL[,"normal.bam"] == 0,]
# somatic loss of heterozygosity has no support in the tumour
#loh_vcf <- vcf[geno(vcf)$QUAL[,"tumour.bam"] == 0,]

# Output BEDPE for use by circos
FR07886867_vcf_example_gr <- breakpointRanges(FR07886867_vcf_example)
FR07886867_vcf_example_gr_bedpe <- data.frame(
  chrom1=seqnames(FR07886867_vcf_example_gr),
  start1=start(FR07886867_vcf_example_gr) - 1,
  end1=end(FR07886867_vcf_example_gr),
  chrom2=seqnames(partner(FR07886867_vcf_example_gr)),
  start2=start(partner(FR07886867_vcf_example_gr)) - 1,
  end2=end(partner(FR07886867_vcf_example_gr)),
  name=names(FR07886867_vcf_example_gr),
  score=FR07886867_vcf_example_gr$QUAL,
  strand1=strand(FR07886867_vcf_example_gr),
  strand2=strand(partner(FR07886867_vcf_example_gr))
)
# Just the lower of the two breakends so we don't output everything twice
FR07886867_vcf_example_gr_bedpe_low <- FR07886867_vcf_example_gr_bedpe[str_detect(FR07886867_vcf_example_gr_bedpe$name, "gridss.+o"),]
write.table(FR07886867_vcf_example_gr_bedpe, "FR07886867_vcf_example_gr_bedpe", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
browseVignettes(package = "StructuralVariantAnnotation")


####GRanges learning
gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               IRanges(1:10, width=10:1))
gr0