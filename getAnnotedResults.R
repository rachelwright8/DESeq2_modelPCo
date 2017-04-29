setwd("~/Documents/aims2013/microbialCommunity/5april_qiimeMin1000/DESeq2/")

gg = read.delim("~/Documents/genomes/amil_moya/amil_iso2gene.tab", header = F)
head(gg)

load("ddsResults_LRT.Rdata")
res = as.data.frame(rr)
head(res)
res$V1 = row.names(res)
annot.res = merge(res, gg, by.y="V1", all.x = T)
head(annot.res)

table(annot.res$padj<0.1) # 123

sig = annot.res[annot.res$padj<0.1 & !is.na(annot.res$padj),]
head(sig)

names(sig)[8] = "annotation"

write.table(sig, "sig_annotated_DESeq2_results_PC1.txt", sep="\t",quote=F, row.names=F)
