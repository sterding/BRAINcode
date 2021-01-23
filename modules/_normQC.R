###########################################
# Rscript to generate the RLE and NUSe plot to detect outlier
# Usage: Rscript $PATH/_normQC.R genes.fpkm.allSamples.uniq.xls genes.fpkm.allSamples.uniq.QC.pdf
# Reference: https://www.genevestigator.com/userdocs/manual/qc.html
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Aug-8
###########################################

args<-commandArgs(TRUE)

if(!require(ape)) install.packages('ape'); 
library(ape)
if(!require(reshape2)) install.packages('reshape2'); 
library(reshape2)


FPKMfile=args[1]  # either filename or stdin
matrix=args[2]
outputfile=args[3]

if(is.na(outputfile)) outputfile="samples.QC.plot.pdf"

# FPKMfile="genes.fpkm.cuffnorm.allSamples.uniq.xls"; outputfile="genes.fpkm.HCILB.uniq.QC.pdf" 
# FPKMfile="genes.fpkm.cufflinks.TCPY.uniq.xls"; outputfile="QC.genes.fpkm.cufflinks.TCPY.uniq.pdf" 

message("loading data...")

fpkm=read.table(file(FPKMfile), header=T, check.names =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; 
if(grepl("cufflinks",FPKMfile)) fpkm=fpkm[,grep("_rep",colnames(fpkm))]; ## using *_rep* to filter sample names. This might be buggy for samples not named in this way.
colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))
colnames(fpkm)=gsub("_0$","",colnames(fpkm))

# filter out stranded, unamplified, _SN_
fpkm=fpkm[, grep("PD_|stranded|unamplified|_SN_", colnames(fpkm), invert = T)]

message(paste("# data dim:", dim(fpkm)))

pdf(outputfile, width=7, height=7)

if(file.exists(matrix)){
  message("generating kmer distance plot...")
  
  N=read.table(matrix, header=F,nrows=1)$V1
  nms=as.character(read.table(matrix, header=F, skip=1, nrows=N)$V1)
  df=data.matrix(read.table(matrix,  fill = TRUE, skip=N+1, col.names = 1:(N-1)))
  upper=rbind(cbind(0,t(df)),0); upper[is.na(upper)]=0
  lower=rbind(0,cbind(df,0)); lower[is.na(lower)]=0
  df=lower+upper; diag(df)=NA; 
  kmer=apply(df,1,function(x) median(x,na.rm=T))
  names(kmer)=nms
  
  hist(kmer, breaks=50, xlab="Median k-mer distance", ylab="Number of samples")
  legend("topright", paste(names(kmer[kmer>=0.0009]), round(kmer[kmer>=0.0009],4), sep=": "), bty='n', cex=0.5)
  
}

message("generating RLE plot...")

# RLE: For each gene and each sample, ratios are calculated between the expression of a gene and the median expression of this gene across all samples of the experiment. For each sample, these relative expression values are displayed as a box plot. Since it is assumed that in most experiments only relatively few genes are differentially expressed, the boxes should be similar in range and be centered close to 0.

#Two effects may characterize arrays with lower quality: 1) the spread is greater than that of other arrays from this experiment, and 2) the box is not centered near 0.

# filter genes with 0 in >90% samples
notAllZero <- (rowMeans(fpkm>0)>0.1)
logfpkm=fpkm[notAllZero,]
logfpkm=log10(logfpkm + 1e-4)  # so row value of 0 will be -2 in the transformed value
rle=logfpkm-apply(logfpkm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
op=par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.2, main="Relative Log Expression", xlab="", ylab="RLE", frame=F)

#The other graphical representation (NUSE) represents normalized standard error (SE) estimates from the PLM fit. The SE estimates are normalized such that for each probe set, the median standard error across all arrays is equal to 1. A box plot of NUSE values is drawn for each array. On the NUSE plot, arrays with lower quality will have boxes that are centered higher and/or have a larger spread than the other good quality arrays from the same experiment. Typically, boxes centered above 1.1 represent arrays that have quality problems which are often detected in several of the other QC measures presented in this chapter.

message("generating clustering plot...")

# clustering
sampleDists = 1 - cor(fpkm, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
#plot(hc, cex=0.7, xlab='', main="Cluster Dendrogram (dis = 1 - Spearman_rank_correlation, linkage = complete)")
celltype=gsub("(.*)_.*_(.*)_.*_.*","\\1_\\2",hc$labels)
# [optional] merge HC_SNDA and ILB_SNDA into HCILB_SNDA 
celltype[celltype=='HC_SNDA']='HCILB_SNDA'; celltype[celltype=='ILB_SNDA']='HCILB_SNDA'

batch=gsub(".*_.*_.*_(\\d+)_.*","\\1",hc$labels)
subject=gsub(".*_(.*)_.*_.*_.*","\\1",hc$labels)

# gsurl='https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/pub?gid=1995457670&output=tsv'
# library(RCurl)
# colorcode=read.delim(textConnection(getURL(gsurl)))

library(googlesheets4)
gsurl="https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/edit#gid=1995457670"
gs4_deauth()
colorcode = as.data.frame(read_sheet(gsurl, sheet = 'color code'))
cols = subset(colorcode,GROUP=="cell type")
celltype.colors=paste0("#",cols$HEX[match(celltype, cols$ITEM)])

hc$labels=subject
tree=as.phylo(hc)

## Update: fix edge.color bug. See https://stackoverflow.com/a/22102420/951718
myLabels <- c('node', sort(unique(batch)))
#myColors <- c("black", gray.colors(length(unique(batch)),start=0))
myColors <- c("black", rainbow(length(unique(batch))))
## match colors and labels (nomatch == node => select idx 1)
## (myLabels are reordered by edge ordering
batchColors <- myColors[match(batch[tree$edge[,2]], myLabels, nomatch=1)]

par(mar=c(1,1,1,1))
plot(tree, type = "unrooted", 
     cex=.3, lab4ut='axial',underscore = T, 
     tip.color=celltype.colors, 
     edge.color= batchColors, 
     main="Clustering of samples based on Spearman correlation")
legend("bottomleft", 
       c("-- cell type --",unique(celltype),"-- batch --",paste("batch",sort(unique(batch)))),
       text.col=c('black',paste0("#",cols$HEX[match(unique(celltype), cols$ITEM)]), myColors), 
       bty='n', cex=.5)

message("generating D-statistic plot...")

# D-statistic
par(op)
D=apply(1-sampleDists, 1, median)
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of D-statistic")
cutoffD=quantile(D, probs = 0.05) # 5% quantitle
if(sum(D<cutoffD)) legend("topleft", paste(names(sort(D[which(D<cutoffD)])), round(sort(D[which(D<cutoffD)]),2)), bty='n', cex=.5)

message("generating gender-match plot...")

# Gender-specific expression
# chrY-unique regions (chrY - pseudoautosomal regions)
# ref: http://useast.ensembl.org/info/genome/genebuild/assembly.html
# grep -w chrY gencode.v19.annotation.bed12 | intersectBed -a - -b chrY.PAR.bed -v -wa | cut -f4 | sed 's/__/\t/g' | cut -f2 | sort -u > chrY-unique.geneID.txt
#chrY=read.table('/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/chrY-unique.geneID.txt')
#chrY=logfpkm[chrY$V1,]
#chrY=apply(chrY, 2, mean)
#chrX=logfpkm[chrX,]
chrX='ENSG00000229807'  # XIST
chrY='ENSG00000129824'  # RPS4Y1

#gsurl='https://docs.google.com/spreadsheets/d/e/2PACX-1vQFQ4aQj0sD9oxIqaZ-cEgo7kWcCmNYGBH9emLw8iNu0f6TTjKE5Lte7IBfoMMy57cLjA4pXE0YlPY2/pub?gid=28&output=tsv'
#covariate=read.delim(textConnection(getURL(gsurl)))
gsurl="https://docs.google.com/spreadsheets/d/1McWI4zAtC1qIS4wiYMXvToTKrDkFXaiAqhfDPDHyNIA/edit#gid=28"
gs4_deauth()
covariate = as.data.frame(read_sheet(gsurl, sheet = 'Subjects_info'))
head(covariate)
sex=covariate$SEX[match(subject, covariate$SOURCE_SUBJECT_ID)]

ind=c(grep(chrX,rownames(logfpkm))[1], grep(chrY,rownames(logfpkm))[1])

d=as.data.frame(t(logfpkm[ind,])); colnames(d)=c("chrX","chrY"); d$SEX=sex;
with(d, plot(chrX, chrY, xlab="Expression of XIST", ylab="Expression of RPS4Y1", col= 'white', bg=ifelse(SEX=="F",'red','blue'), pch=21, bty="n", main="Gender-specific expression"))
#if(nrow(subset(d, chrX>0 & chrX<1.5 & chrY<0.2))>0) text(subset(d, chrX>0 & chrX<1.5 & chrY<0.2), rownames(subset(d, chrX>0 & chrX<1.5 & chrY<0.2)),pos=2, cex=0.5)
mean2sd_chrX=mean(d$chrX[sex=="F"])-sd(d$chrX[sex=="F"])
mean2sd_chrY=mean(d$chrY[sex=="M"])-sd(d$chrY[sex=="M"])
if(nrow(subset(d, (SEX=="F" & chrX<mean2sd_chrX) | (SEX=="M" & chrY<mean2sd_chrY)))>0) text(subset(d, (SEX=="F" & chrX<mean2sd_chrX) | (SEX=="M" & chrY<mean2sd_chrY)), rownames(subset(d, (SEX=="F" & chrX<mean2sd_chrX) | (SEX=="M" & chrY<mean2sd_chrY))),pos=2, cex=0.5)
legend('bottomleft',pch=21,c("Female","Male"), col='white',pt.bg=c("red","blue"), bty='n', cex=1)

dev.off()

message("QC done")
