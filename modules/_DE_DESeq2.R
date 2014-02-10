###########################################
# R script for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 1/27/2014
# version: 1.0
# Note: 
###########################################

## TODO: check if any internal filter for low read count genes? or gene only expressed in few samples?


input_dir="~/neurogen/rnaseq_PD/run_output/"

## read covariance table
## via remote URL (Tips: http://blog.revolutionanalytics.com/2009/09/how-to-use-a-google-spreadsheet-as-data-in-r.html)
## may not work in cluster node
#covarianceTableURL="https://docs.google.com/spreadsheet/pub?key=0Aumm3V3g3dF7dEFnZ2pPQjlheXlZand6YWUxeF9PMUE&single=true&gid=5&output=csv"  # url for covariance data
#require(RCurl)
#covarianceTable=read.csv(textConnection(getURL(covarianceTableURL)))

# via local file
covarianceTableURL="~/neurogen/rnaseq_PD/results/DE_DESeq2/covariances.tab"
covarianceTable=read.table(covarianceTableURL, header=T)

# TODO: install packages if not available
# source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")
library("DESeq2")

###########################################
# step1: load co-variances tabel
###########################################

#sampleFiles <- list.files(path=input_dir, pattern="^(ILB|PD)")
sampleFiles <- list.files(path=input_dir, pattern="^(PD|HC|ILB)")
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = paste(sampleFiles, "hgseqcount.by.gene.tab", sep="/"))  # tocheck


# subset and re-order
covarianceTable$sampleName=as.character(covarianceTable$sampleName)
covarianceTable=covarianceTable[match(intersect(sampleFiles, covarianceTable$sampleName), covarianceTable$sampleName),]
rownames(covarianceTable)=covarianceTable$sampleName;
covarianceTable=covarianceTable[,-1]
covarianceTable$condition=factor(as.character(covarianceTable$condition))
sampleTable=cbind(sampleTable, covarianceTable)

###########################################
# step2: load HTSeq output
###########################################
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = input_dir,
                                       design= ~ condition + RIN + PMI + batch + age + sex)

#colData(dds)$condition <- factor(colData(dds)$condition, levels=c("ILB", "PD"))
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("HC","ILB", "PD"))

###########################################
# step3: diagnosis of the data [optional]
###########################################

# note: instead of using rank, use mean directly
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

pdf("diagnosis.pdf", width=15, height=5)
par(mfrow=c(1,3))
plot(rank(rowMeans(counts(dds))), genefilter::rowVars(counts(dds)),
main="no transformation", ylim=c(0,3000))
plot(rank(rowMeans(counts(dds))), genefilter::rowVars(log2(counts(dds)+1)),
main="log2(x+1) transform")
plot(rank(rowMeans(assay(vsd))), genefilter::rowVars(assay(vsd)),
main="VST")

dev.off()

cds=dds
design(cds) <- ~ 1
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
vsd <- getVarianceStabilizedData( cds )

## save the variance-stabilized data
write.table(vsd, "htseqcount.vst.allsamples.xls", sep="\t", quote = F, col.names = NA, row.names = TRUE)

pdf("hist.allsample.pdf")
#sapply(1:ncol(vsd), function(x) {hist(vsd[,x], main=colnames(vsd)[x])})
t=sapply(1:ncol(vsd), function(x) {hist(vsd[,x], xlab="variance stabilized count", main=colnames(vsd)[x], breaks=100)})
dev.off()

pdf("hist2.allsample.pdf")
#sapply(1:ncol(vsd), function(x) {hist(vsd[,x], main=colnames(vsd)[x])})
t=sapply(1:ncol(vsd), function(x) {xx=vsd[,x]; hist(xx[xx>4], xlab="variance stabilized count", main=colnames(vsd)[x], breaks=100)})
dev.off()

## scatterplot of all sample pairs
## put histograms on the diagonal
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
panel.smooth2 <- function (x, y) 
{
    points(x, y, col = "#0000ff20", pch = 20)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) abline(lm(y~x),  col = 'red')
}

#pairs(vsd[,1:3], lower.panel=panel.smooth2, upper.panel=panel.cor)

pdf("scatterplot.pdf", height=30, width=30)
library(gclus)
dta <- vsd
dta.r <- abs(cor(dta)) # get correlations
dta.col <- dmat.color(dta.r) # get colors
# reorder variables so those with highest correlation
# are closest to the diagonal
dta.o <- order.single(dta.r) 
#cpairs(dta, dta.o, panel.colors=dta.col, gap=.5, main="Variables Ordered and Colored by Correlation" )
cpairs(dta, dta.o, lower.panel=panel.smooth2, upper.panel=panel.cor, panel.colors=dta.col, gap=.5, main="Variables Ordered and Colored by Correlation" )
dev.off()

## Hierarchical clustering of samples, using correlations between variables "as distance"

pdf("clustering.tree.pdf", width=15, height=5)
plot(hclust(as.dist((1 - cor(vsd))),method = "single"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=1-correlation, linkage=single)") 
plot(hclust(as.dist((1 - cor(vsd))),method = "complete"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=1-correlation, linkage=complete)")
plot(hclust(as.dist((1 - cor(vsd))),method = "ave"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=1-correlation, linkage=average)")
plot(hclust(as.dist((1 - cor(vsd))),method = "ward"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=1-correlation, linkage=ward)")
plot(hclust(dist(t(vsd)),method = "single"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=euclidean, linkage=single)")
plot(hclust(dist(t(vsd)),method = "complete"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=euclidean, linkage=complete)")
plot(hclust(dist(t(vsd)),method = "average"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=euclidean, linkage=average)")
plot(hclust(dist(t(vsd)),method = "ward"), cex=0.7, xlab='', main="Cluster Dendrogram (dis=euclidean, linkage=ward)")
dev.off()

# heatmap
library("gplots")
heatmap.2(vsd[1:100,1:10],
dendrogram="column",
Rowv=NULL,
Colv=TRUE,
distfun = function(x) as.dist(1 - cor(x)),
hclustfun = function(x) hclust(x,method = 'ave'))

heatmap.2(mat,
dendrogram="column",
Rowv=NULL,
Colv=TRUE,
distfun = function(x) as.dist(1 - cor(x)),
hclustfun = function(x) hclust(x,method = 'ave'))

###########################################
# step3: call differential expressed genes
###########################################

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
head(res)

###########################################
# step4:  Exploring and exporting results
###########################################
plotMA(dds, main="DESeq2", ylim=c(-2,2))
write.csv(as.data.frame(res), file="condition_treated_results.csv")


