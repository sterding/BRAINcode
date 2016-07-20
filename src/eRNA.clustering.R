##############################################
## R script to cluster eRNA 
## Author: Xianjun Dong
## Usage: Rscript $HOME/neurogen/pipeline/RNAseq/src/eRNA.clustering.R `ls eRNA.f*.txt`
## Date:
## Version: 0.0
##############################################

setwd("~/eRNAseq")

# ================================
## sample clustering 
# ================================

library(gsheet) # install.packages('gsheet') # covaraice table URL
colorcode = gsheet2tbl('docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118')

## TODO
## 1. clustering of HITNEs using HITNE RPM, DNase, H3K4me1, H3k27ac from Roadmap, and CAGE from FANTOM, in all 5 cell types --> similar as the Fig 3a,b,c of Nature FANTOM5
## 2. aggregation plot for the same marks
## 3. clustering dendrogram accross samples
## 4. heatmap of cor() for all HiTNEs across 90 SNDA samples
df=read.table("~/eRNAseq/eRNA.90samples.meanRPM.xls", header=T)
rownames(df) = df[,1]; df=df[,-1]

## take the one overlapped wit DNase peaks (N=5669)
core=read.table("eRNA.wDNase.bed")

df=df[rownames(df) %in% core$V4,]
df=scale(df)
require(gplots);
heatmap.2(1-cor((df)), trace="none", density="none", scale='none', cexRow=0.5, labCol ="", margins=c(3,8), symm=T, lhei=c(1,5))
#heatmap.2(d, trace="none", density="none", scale='none', cexRow=0.5, labCol ="", labRow ="", margins=c(3,8), symm=T, lhei=c(1,5))

x=cor(t(df), method='spearman')
# remove row/col which max(abs)<0.7
x2=x; x2[lower.tri(x2, T)]=0
INDEX=apply(abs(x2),1,max, na.rm=T)>0.7
x=x[INDEX,INDEX]

pdf("eRNA.wDNase.cor.clustering.pdf")
heatmap.2(x, breaks=seq(-1,1,0.01),col=colorRampPalette(c('darkblue','white','red'))(200), trace="none", density="none", symm=F,symkey=F,symbreaks=T, scale='none', cexRow=0.5, labCol ="", margins=c(3,8), lhei=c(1,5))
dev.off()

xmin=-1; xmax=1;
x[x<xmin]=xmin; x[x>xmax]=xmax;
x=x2
x[lower.tri(x)]=0

#collist<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
collist<-c('darkblue','white','red')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]

pdf("eRNA.wDNase.cor.image.pdf")
par(mar=c(0,2,0,2), oma=c(3,3,3,3))
layout(matrix(1:2,nrow = 2), widths = 8, heights = c(1,8))
require('zoo'); # install.packages('zoo')
r=rollapply(apply(x2,2,mean),width =5,by = 1,FUN = mean, align = "left")
plot(r, type='h', col='gray', xlab="",ylab="",xaxs = 'i', xaxt="n",yaxt="n", frame.plot=F)
image(as.matrix(x), col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=F)
dev.off()


df=read.table("~/eRNAseq/eRNA.140samples.meanRPM.xls", header=T)
rownames(df) = df[,1]; df=df[,-1]
df=df[,grep("[HC|ILB]_.*_.*", colnames(df))]

cols = subset(colorcode,GROUP=="cell type")
CL=paste0("#",cols$HEX[match(gsub("(.*)_.*_(.*)_[1-6]_rep.*","\\1_\\2",colnames(df)), cols$ITEM)])
# scross samples
require(gplots);
heatmap.2(1-cor(df, method='spearman'), trace="none", density="none", sepwidth=c(0,0), sepcolor=NA, scale='none', cexRow=0.4, labCol ="", margins=c(3,8), symm=T, lhei=c(1,5), RowSideColors=CL)

# scross HITNEs
df0=df[grep("chr18|chr19|chr20|chr21|chr22|chrX|chrY", rownames(df)),]
# sort by chr-->start
d=data.frame(do.call(rbind, strsplit(rownames(df0),split="_")), stringsAsFactors =F)
d$X2=as.numeric(d$X2); d$X3=as.numeric(d$X3)
df0=df0[with(d,order(X1,X2)),]
rownames(df0)
#x=cor(t(df0), method='spearman')
x=cor(t(log10(0.001+df0)), method='pearson')

xmin=-1; xmax=1; x[x<xmin]=xmin; x[x>xmax]=xmax;
collist<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
ColorRamp<-colorRampPalette(collist)(10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]

par(mar=c(0,2,0,2), oma=c(3,3,3,3))
layout(matrix(1:3,nrow = 3), widths = 8, heights = c(1,8,1))
r=rollapply(apply(x,2,mean),width =5,by = 1,FUN = mean, align = "left")
firstm=match(levels(factor(gsub("(.*)_.*_.*","\\1",names(r)))), gsub("(.*)_.*_.*","\\1",names(r)))
lastm=rev(length(r)+1-match(rev(levels(factor(gsub("(.*)_.*_.*","\\1",names(r))))), rev(gsub("(.*)_.*_.*","\\1",names(r)))))
mid=(firstm+lastm)/2
plot(r, type='h', col=1+match(gsub("(.*)_.*_.*","\\1",names(r)),levels(factor(gsub("(.*)_.*_.*","\\1",names(r))))), xlab="",ylab="",xaxs = 'i', xaxt="n",yaxt="n", frame.plot=F)
mtext(levels(factor(gsub("(.*)_.*_.*","\\1",names(r)))), 1, line=-2, at=mid)
image(as.matrix(x), col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=F)
r=rollapply(apply(df0,1,mean),width =5,by = 1,FUN = mean, align = "left")
plot(-r, type='h', col=1+match(gsub("(.*)_.*_.*","\\1",names(r)),levels(factor(gsub("(.*)_.*_.*","\\1",names(r))))), xlab="",ylab="",xaxs = 'i', xaxt="n",yaxt="n", frame.plot=F)


## 5. heatmap of cor() for all HiTNEs across different cell types
## 6. heatmap of cor() for all HiTNEs x GENCODE genes

quit('no')


expressionFile=args[1]  # e.g. eRNAfinal.allsamples.RPKM.tab

filename=sub("(.*)\\.[^.]*", "\\1",expressionFile) # e.g. eRNAfinal.allsamples.RPKM

df0=read.table(expressionFile, header=T)
rownames(df0)=df0[,1]; df0=df0[,-1]
df=log10(df0+1)  # log transformed

# filter out samples
hist(rowMeans(df), breaks=100)
df=df[,colMeans(df==0)<0.5]

# cluster the samples/columns
# ----------------------
hc=hclust(as.dist((1 - cor(df))))  # cluster of samples (dis=1-correlation, linkage=complete)

pdf(paste(filename, "colclusters","pdf", sep="."))
plot(hc, cex=0.5, xlab='')
library(ape)
plot(as.phylo(hc), type = "unrooted", lab4ut='axial', cex=.5)
dev.off()

# re-order the columns
df=df[,hc$order]

# cluster the genes/rows
# ----------------------
set.seed(1)
# decide of the number of k-mean clusters to build
# TODO: use clValid or other functions to better define this number
maxclust <- 10
y=numeric(maxclust);
for(i in 1:maxclust) {
    cat(i)
    cl=kmeans(df, i)
    y[i]=100*(1-cl$tot.withinss/cl$totss);
}
pdf(paste(filename, "kmeans.cluster","pdf", sep="."))
plot(y, xlab="Number of clusters", ylab="% of explained variance", type='b')
dev.off()
maxclust <- 5

cl<- kmeans(df,maxclust)

# re-order the rows
o<- order(cl$cluster, -rowMeans(df)) # order rows by the kmean cluster, and then the rowMean
df<- df[o,]

library(pheatmap) # IF NOT, install.packages('pheatmap')
library("RColorBrewer")
# create heatmap

# Generate column annotations
annotation=data.frame(do.call(rbind, strsplit(colnames(df),"_")))
rownames(annotation) = colnames(df)
colnames(annotation) = c("Diagnosis","SubjectID","Cell_Type","Batch")
annotation=annotation[,c("Diagnosis","Cell_Type","Batch")]

# Specify colors for annotation
cHC=rgb(44,162,95,maxColorValue =255); cILB=rgb(254,178,76,maxColorValue =255); cPD=rgb(255,0,0,maxColorValue =255)
cDiagnosis = c(cHC,cILB,cPD); names(cDiagnosis) = c("HC", "ILB", "PD")
cCell_Type = c(rgb(44,162,95,maxColorValue =255), rgb(153,216,201,maxColorValue =255),rgb(229,245,249,maxColorValue =255)); names(cCell_Type) = c("SNDA", "MCPY", "TCPY")
cBatch = colorRampPalette(c("lightgray", "black"))(length(levels(annotation$Batch))); names(cBatch)=levels(annotation$Batch)
ann_colors = list(Diagnosis = cDiagnosis, Cell_Type = cCell_Type, Batch = cBatch)

range(df)
hist(rowMeans(df), breaks=100)
df[df>4]=4
bk=unique(c(seq(0,2, length=100),seq(2,4,length=100)))
hmcols<- colorRampPalette(c('blue',"white","red"))(length(bk)-1)  

hm.parameters <- list(df,
  breaks=bk,
  color = hmcols,
  cellwidth = 8, cellheight = 8,
  border_color=NA,
  fontsize=7,
  scale = "none",
  kmeans_k = 8,
  show_rownames = T, show_colnames = T,
  main = "Heatmap of eRNA clustering",
  clustering_method = "complete",
  cluster_rows = F, 
  clustering_distance_rows = "euclidean",
  cluster_cols = T,
  clustering_distance_cols = dist(t(df)), #as.dist(1-cor(df)),
  treeheight_col = 30,
  annotation = annotation,
  annotation_colors = ann_colors
)

## To draw the heat map on screen 
#do.call("pheatmap", hm.parameters)

# To draw to file 
do.call("pheatmap", c(hm.parameters, filename=paste(filename, "hm.rowclusters","pdf", sep=".")))

hm.parameters$kmeans_k=NA; hm.parameters$cellwidth=NA; hm.parameters$cellheight=NA; hm.parameters$show_rownames = F;
# To draw to file 
do.call("pheatmap", c(hm.parameters, filename=paste(filename, "hm.allrows","png", sep=".")))