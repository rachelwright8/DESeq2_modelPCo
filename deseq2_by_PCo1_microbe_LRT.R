setwd("~/Documents/aims2013/microbialCommunity/5april_qiimeMin1000/DESeq2/")
library(DESeq2)

# load microbial population data -----------------------------------------------
data16s0 = read.csv("~/Documents/aims2013/microbialCommunity/5april_qiimeMin1000/pcoa_vegan/weighted_forPCoA.csv", header=F)
head(data16s0)
row.names(data16s0) <- data16s0$V1
data16s0$V1 <- NULL
head(data16s0)
# round to 3 decimal places
data16s = round(data16s0,3)
head(data16s)
summary(data16s)

# subset for first three principal components
data16s = data16s[c(1:3)]
head(data16s)
names(data16s) = c("PC1", "PC2", "PC3")
head(data16s)
# the "1" from this sample name got lost. replace so the sample is matched with other datasets.
row.names(data16s)[15] = "L71"

PCo1 = data16s[c(1)]
head(PCo1)
summary(PCo1)
hist(PCo1,3)
# low survival corals tended to have higher PCo1 values
# gene expression will be modeled as LFC / unit PCo1 value
# upregulated = highly expressed in samples on the right side of the PCo1 axis (tend to be low survival)
# downregulated = highly expressed in samples on the left side of hte PCo1 axis (tend to be high survival)

# Load data ---------------------------------------------------------------
countdata = read.table("~/Documents/aims2013/RNAseq/Good_allcounts_aims2013_april2015.txt",header=TRUE,row.names = 1) 
head(countdata) 
length(countdata[,1])
# 44687 genes mapped

names(countdata)
names(countdata) = sub("*.fq.trim.sam.counts","",names(countdata))
names(countdata) = gsub("[.]","",names(countdata)) 

# subset for samples in microbial community dataset
# any samples missing from the countdata?
data16s[!row.names(data16s) %in% names(countdata),]
# no!

counts = countdata[,names(countdata) %in% row.names(data16s)]
head(counts)
mean(colSums(counts)) #405382.8 counts / sample

# Set base mean minimum ---------------------------------------------------
means = apply(counts,1,mean)
table(means>3)
# FALSE  TRUE 
# 29551 15136  

means3 = names(means[means>3])
head(means3)
length(means3)
#15136

countDataFilt = counts[row.names(counts) %in% means3,]
head(countDataFilt)

totalCountsFilt = colSums(countDataFilt)
totalCountsFilt

mean(totalCountsFilt)#385544.5
max(totalCountsFilt) #750002
mean(totalCountsFilt) #385544.5

# add sample data ------------------------------------------------------------
coldata0 = read.delim("~/Documents/aims2013/microbialCommunity/5april_qiimeMin1000/all.qiime.map.txt")
coldata0 = coldata0[c(1,5:9)]
head(coldata0)
row.names(coldata0) = coldata0$X.SampleID
coldata0$X.SampleID = NULL
head(coldata0)
row.names(coldata0)[11] = "L71" # fix this sample name again
coldata = coldata0[match(names(counts), row.names(coldata0)),]
head(coldata)

# any samples missing from the sample data?
data16s[!row.names(data16s) %in% row.names(coldata0),]
# no!

# match sample order
countDataFilt = countDataFilt[names(countDataFilt) %in% row.names(coldata)]
head(countDataFilt)
# double-check sample order (matters for deseq2!!)
identical(names(countDataFilt),row.names(coldata))
#TRUE

# add PCo1 to coldata
coldataPC = merge(coldata, PCo1, by=0)
head(coldataPC)

# Construct data object ---------------------------------------------------

dds1 <- DESeqDataSetFromMatrix(
  countData = countDataFilt,
  colData = coldataPC,
  design = ~ PC1 + bacYN + surv)

# save(coldataPC, countDataFilt, dds1, file="dds_LRT.Rdata")

# Run LRT model -----------------------------------------------------------
deds1 = DESeq(dds1,test="LRT",reduced=~ bacYN + surv)

# Extract results ---------------------------------------------------------
rr = results(deds1, name = "PC1")
head(rr)
# log2 fold change (MLE): PC1 
# LRT p-value: '~ PC1 + bacYN + surv' vs '~ bacYN + surv' 

table(rr$log2FoldChange>0)
# FALSE  TRUE 
# 77694  7442 

# correct the sign of the statistic
rr$stat[rr$log2FoldChange<0] = (-1)*rr$stat[rr$log2FoldChange<0]

# how many genes are significantly downregulated (increasing expresion from right to left?)
table(rr$padj<0.1 & rr$log2FoldChange<0)
# FALSE  TRUE 
# 14737    72 

# how many genes are significantly downregulated (increasing expresion from left to right?)
table(rr$padj<0.1 & rr$log2FoldChange>0)
# FALSE  TRUE 
# 14825    51 

# 123 genes DE (FDR = 0.1)

# Get RLD (this will take awhile) --------------------------------------------------------
rld = rlogTransformation(deds1)

# Make rlogdata and pvals table
valsPC1 = cbind(rr$pvalue, rr$padj)
colnames(valsPC1)=c("pval.p1", "padj.p1")
rldpvals=cbind(assay(rld),valsPC1)
head(rldpvals)
dim(rldpvals)
#15136    33
# write.csv(rldpvals, "rld_pvals_by_PC1_LRT.csv", quote=F)

# save(coldataPC, countDataFilt, dds1, deds1, rr, rld, file="ddsResults_LRT.Rdata")

# heatmap ------------------
source("~/Documents/scripts/uniHeatmap.R") # for heatmap

gg = read.delim("~/Documents/genomes/amil_moya/amil_iso2gene.tab", header = F)
head(gg)

exp = read.csv("rld_pvals_by_PC1_LRT.csv")
head(exp)
row.names(exp) = exp$X
exp$X = NULL
names(exp)
# get rid of columns with p values
rldd = exp[-c(32:33)]
head(rldd)

quartz()
nums = uniHeatmap(vsd = rldd, gene.names = gg, sort = coldataPC$PC1, metric = exp$padj.p1, cutoff = 0.05, pdf = F,cex = 0.9, pval = 0.05)
nums 

# Write results for GO/KOG analysis -------------------------------------------

# by STAT
stat=data.frame(cbind("gene"=row.names(rr),"stat"=rr$stat))
write.csv(stat,quote=F,row.names=F,file="GO_PC1_STAT.csv")

# by LFC
lfc=data.frame(cbind("gene"=row.names(rr),"lfc"=rr$log2FoldChange))
write.csv(lfc,quote=F,row.names=F,file="GO_PC1_LFC.csv")

