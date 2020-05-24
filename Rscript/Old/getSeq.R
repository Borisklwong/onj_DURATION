#load hg19 reference genome from USSC
library(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)
# index of genome
Hsapiens
# show all genome names
seqnames(Hsapiens)
#get the sequence chr17:30276680-30276880
getSeq(Hsapiens, "chr17", 30276680, 30276880)



gridss105bb_4820_1 = getSeq(Hsapiens, "chr17", 30276680, 30276880)
gridss105bb_4820_2 = getSeq(Hsapiens, "chr16", 72888843, 72889043)

gridss105bb_4820_2
