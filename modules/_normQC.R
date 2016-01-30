###########################################
# Rscript to generate the RLE and NUSe plot to detect outlier
# Usage: Rscript $PATH/_normQC.R genes.fpkm.allSamples.uniq.xls genes.fpkm.allSamples.uniq.QC.pdf
# Reference: https://www.genevestigator.com/userdocs/manual/qc.html
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Aug-8
###########################################

args<-commandArgs(TRUE)

require(ape)
require(reshape2)

FPKMfile=args[1]  # either filename or stdin
outputfile=args[2]

# FPKMfile="genes.fpkm.HCILB.uniq.xls"; outputfile="genes.fpkm.HCILB.uniq.QC.pdf" 

message("loading data...")

fpkm=read.table(file(FPKMfile), header=T, check.names =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; fpkm=fpkm[,grep("FPKM",colnames(fpkm))]; colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))


message("generating RLE plot...")

# RLE: For each gene and each sample, ratios are calculated between the expression of a gene and the median expression of this gene across all samples of the experiment. For each sample, these relative expression values are displayed as a box plot. Since it is assumed that in most experiments only relatively few genes are differentially expressed, the boxes should be similar in range and be centered close to 0.

#Two effects may characterize arrays with lower quality: 1) the spread is greater than that of other arrays from this experiment, and 2) the box is not centered near 0.

pdf(outputfile, width=7, height=7)
par(mfrow=c(2,2))
# filter genes with 0 in >90% samples
notAllZero <- (rowMeans(fpkm>0)>0.1)
logfpkm=fpkm[notAllZero,]
logfpkm=log10(logfpkm + 1e-4)  # so row value of 0 will be -2 in the transformed value
rle=logfpkm-apply(logfpkm, 1, median) # change / to - so that we got log(fold-change) which centered on 0 on the RLE plot.
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
op=par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="Relative Log Expression", xlab="", ylab="RLE", frame=F)

#The other graphical representation (NUSE) represents normalized standard error (SE) estimates from the PLM fit. The SE estimates are normalized such that for each probe set, the median standard error across all arrays is equal to 1. A box plot of NUSE values is drawn for each array. On the NUSE plot, arrays with lower quality will have boxes that are centered higher and/or have a larger spread than the other good quality arrays from the same experiment. Typically, boxes centered above 1.1 represent arrays that have quality problems which are often detected in several of the other QC measures presented in this chapter.

message("generating clustering plot...")

# clustering
sampleDists = 1 - cor(fpkm, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
#plot(hc, cex=0.7, xlab='', main="Cluster Dendrogram (dis = 1 - Spearman_rank_correlation, linkage = complete)")
celltype=gsub(".*_.*_(.*)_.*_.*","\\1",hc$labels)
batch=gsub(".*_.*_.*_(.*)_.*","\\1",hc$labels)
hc$labels=gsub(".*_(.*)_.*_.*_.*","\\1",hc$labels)
par(op)
plot(as.phylo(hc),type = "unrooted", cex=.5, lab4ut='axial',underscore = T, tip.color=ifelse(celltype=="SNDA",rgb(44,162,95,maxColorValue =255),ifelse(celltype=="TCPY",rgb(158,188,218,maxColorValue =255),rgb(153,216,201,maxColorValue =255))), edge.color= gray.colors(5,start=0)[6-as.numeric(batch)], main="Clustering of samples based on Spearman correlation")
legend("bottomleft", c("-- cell type --","SNDA","MCPY","TCPY","-- batch --",paste("batch",1:5)),text.col=c('black', rgb(44,162,95,maxColorValue =255), rgb(153,216,201,maxColorValue =255), rgb(158,188,218,maxColorValue =255),'black',gray.colors(5,start=0)[5:1]), bty='n')

message("generating D-statistic plot...")

# D-statistic
D=apply(1-sampleDists, 1, median)
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of D-statistic")
legend("topleft", paste(names(sort(D[which(D<0.7)])), round(sort(D[which(D<0.7)]),2)), bty='n')

message("generating gender-match plot...")

# Gender-specific expression
# chrY-unique regions (chrY - pseudoautosomal regions)
# ref: http://useast.ensembl.org/info/genome/genebuild/assembly.html
# grep -w chrY gencode.v19.annotation.bed12 | intersectBed -a - -b chrY.PAR.bed -v -wa | cut -f4 | sed 's/__/\t/g' | cut -f2 | sort -u > chrY-unique.geneID.txt
#chrY=read.table('/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/chrY-unique.geneID.txt')
#chrY=logfpkm[chrY$V1,]
#chrY=apply(chrY, 2, mean)
#chrX=logfpkm[chrX,]
chrX='ENSG00000229807.5'  # XIST
chrY='ENSG00000129824.11'  # RPS4Y1
covariate=read.table("~/neurogen/rnaseq_PD/rawfiles/covariances.tab", header=T)
rownames(covariate)=covariate[,1]
d=as.data.frame(t(logfpkm[c(chrX,chrY),])); colnames(d)=c("chrX","chrY")
plot(d, xlab="Expression of XIST", ylab="Expression of RPS4Y1", col= ifelse(covariate[colnames(fpkm),'sex']=="F",'red','blue'), pch=15, bty="n", main="Gender-specific expression")
text(subset(d, chrX>0 & chrX<1.5 & chrY<0.2), rownames(subset(d, chrX>0 & chrX<1.5 & chrY<0.2)),pos=2)
legend('bottomleft',pch=15,c("Female","Male"), col=c("red","blue"), bty='n', cex=1.5)

dev.off()

message("QC done")
