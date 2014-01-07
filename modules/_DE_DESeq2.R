###########################################
# R script for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 1.0
# Note: 
###########################################

input_dir="~/neurogen/rnaseq_PD/run_output/"
covarianceTableURL="~/neurogen/rnaseq_PD/results/DE_DESeq2/covariances.tab"  # url for covariance data



# TODO: install packages if not available
# source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")
library("DESeq2")

# step1: load co-variances tabel

sampleFiles <- list.files(path=input_dir, pattern="^(PD|HC)")
#sampleFiles <- list.files(path=input_dir, pattern="^(PD|HC|ILB)")
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = paste(sampleFiles, "hgseqcount.by.gene.tab", sep="/"))  # tocheck

# read covariance table
covarianceTable=read.table(covarianceTableURL, header=T)
# subset and re-order
covarianceTable$sampleName=as.character(covarianceTable$sampleName)
covarianceTable=covarianceTable[match(intersect(sampleFiles, covarianceTable$sampleName), covarianceTable$sampleName),]
rownames(covarianceTable)=covarianceTable$sampleName;
covarianceTable=covarianceTable[,-1]
covarianceTable$condition=factor(as.character(covarianceTable$condition))
sampleTable=cbind(sampleTable, covarianceTable)

# step2: load HTSeq output
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = input_dir,
                                       design= ~ condition + RIN + PMI + batch + age + sex)

colData(dds)$condition <- factor(colData(dds)$condition, levels=c("HC", "PD"))
#colData(dds)$condition <- factor(colData(dds)$condition, levels=c("HC","ILB", "PD"))

# step3: call differential expressed genes

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
head(res)

# step4:  Exploring and exporting results
plotMA(dds, main="DESeq2", ylim=c(-2,2))
write.csv(as.data.frame(res), file="condition_treated_results.csv")

