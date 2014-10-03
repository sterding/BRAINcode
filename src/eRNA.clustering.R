##############################################
## R script to cluster eRNA 
## Author: Xianjun Dong
## Usage:
## Date:
## Version: 0.0
##############################################

args<-commandArgs(TRUE)

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