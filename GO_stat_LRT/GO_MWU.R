setwd("~/Documents/aims2013/microbialCommunity/5april_qiimeMin1000/DESeq2/GO_stat_LRT/")

# Edit these to match your data file names ----
goAnnotations="amil_defog_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")

input="GO_PC1_STAT.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).

# Run here -------------------
goDivision="MF" # either MF, or BP, or CC
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", 
	largest=0.1,  
	smallest=5,   
	clusterCutHeight=0.25,
  Alternative="t"
)

# Do not continue if the printout shows that no GO terms pass 10% FDR.

# PC1
# CC 4
# BP 15
# MF 14

# Plotting results --------
quartz()
gomwuPlot(input,goAnnotations,goDivision,
	absValue=2,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.5 if you are doing Fisher's exact test for standard GO enrichment.
	level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
	level2=0.05, # FDR cutoff to print in regular (not italic) font.
	level3=0.01, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.3 # height of the hierarchical clustering tree
)

# get results ---------------
sig.cc = read.table("MWU_CC_GO_PC1_STAT.csv", header=TRUE, quote="\"")
sig.cc[sig.cc$p.adj<0.1,]
# all positive