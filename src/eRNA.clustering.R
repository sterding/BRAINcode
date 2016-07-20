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


# ================================
## HiTNE clustering
# ================================

features = read.table("~/eRNAseq/eRNA.characterize.xls", header=T)

# ---------------------------------------------
# clustering only based on if overlapping with nTFBS, P300, CAGE, chromHMM, DNase, HCNE
# ---------------------------------------------

df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
#df$f20.nHostgene=ifelse(df$f20.nHostgene>=10,1,0)
df[df>0]=1;
library(plyr); 
df0=count(df, vars = colnames(df))
mat_data=df0[,1:6]
#colnames(mat_data)=c("TFBS_hotspot","P300","CAGE_enhancer","")

# do the clustering
require(ade4); require(cluster)
#row_distance = dist.binary(mat_data, method = 10)
row_distance = daisy(mat_data, metric = "gower", type =list(asymm=c(1:6)), weights = c(1,1,1,1,5,1))
row_cluster = hclust(row_distance, method = "single")
pdf("~/Dropbox/PDBrainMap/figures/eRNA/eRNA.clustering.tree.pdf")
plot(row_cluster)
dev.off()

col_distance = as.dist((1 - cor(mat_data))/2)
col_cluster = hclust(col_distance, method = "complete")

gr.row <- cutree(row_cluster, 7)
gr.col <- cutree(col_cluster, 4)

# below code are partially borrowed from http://sebastianraschka.com/Articles/heatmaps_in_r.html
# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
# gr.row <- cutree(row_cluster, 7)
# gr.col <- cutree(col_cluster, 4)
# gr.row[row_cluster$order]=c(1,rep(2,6),rep(3,15),rep(4,13),rep(5,10),rep(6,3))
# col1 <- brewer.pal(6, "Set1")
# col2 <- brewer.pal(4, "Pastel1")

# reorder to have DNase positioned at the left and HiTNEs with none support at the top
require(ape)
row_cluster <- rev(reorder(as.dendrogram(row_cluster), 48:1, min))
col_cluster <- rev(reorder(as.dendrogram(col_cluster), 1:6))
plot(row_cluster)

require(gplots);
require(RColorBrewer);

#col1 <- brewer.pal(5, "Set1"); gr.row[order.dendrogram(row_cluster)]=c(rep(1,22),rep(2,10),rep(3,10),rep(4,5), 5)
col1=c("blue","purple","red"); gr.row[order.dendrogram(row_cluster)]=c(1,rep(2,25), rep(3,22))
col2 <- brewer.pal(4, "Pastel1"); gr.col[order.dendrogram(col_cluster)]=c(1,1,2,2,3,4)

pdf("~/Dropbox/PDBrainMap/figures/eRNA/eRNA.clustering.3class.pdf", width=2.5, height=7)
pdf("eRNA.clustering.3class.pdf", width=2.5, height=7)
heatmap.2(as.matrix(mat_data),
          margins=c(10,6), keysize=0.1, sepcolor=NA,sepwidth=c(0,0),
          lwid=c(1.5,0.25,4), lhei=c(0.1,0.01,0.5),
          #lmat=rbind(c(5,0,4),c(3,1,2)),
          lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3)),
          #cellnote = mat_data,  # same data set for cell labels
          main = NA, # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=colorRampPalette(c("gray", "black"))(n = 2),       # use on color palette defined earlier
          cexCol=1,
          cexRow=1,
          labRow = df0$freq, #[rownames(mat_data)],
          labCol = sub("f[0-9]+.(.*)","\\1",colnames(mat_data)),
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          Colv = as.dendrogram(col_cluster),
          Rowv = as.dendrogram(row_cluster),
          RowSideColors=col1[gr.row],
          ColSideColors=col2[gr.col],
          key.title=NA, key.ylab=NA, key.xlab = NA, key.par=list(mgp=c(3,0.3,0))
)

dev.off()

# ---------------------------------------------
# clustering based on mixed features
# ---------------------------------------------

# mixed features
df=subset(features, select=c(f01.dis2TSS, f02.RPKM, f03.RPM, f05.CpG, f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f11.DNase, f12.DNaseENCODE, f12.DNaseROADMAP,f13.phyloP, f15.HCNE, f16.GWAS, f18.eSNP, f20.nHostgene, f21.lenHostgene, f22.lenHostgeneMetaintron))
df=df[grep("chr7", rownames(df)),]
## divide into continous and discrete variables
# continous
df1=subset(df, select=c(f02.RPKM, f03.RPM, f05.CpG, f06.TFBS, f11.DNase, f20.nHostgene)) #, f21.lenHostgene, f22.lenHostgeneMetaintron))
df1=log10(df1+0.5)
df1=cbind(df1, f13.phyloP = df$f13.phyloP)
# discrete
df2=subset(df, select=c(f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
df2[df2>0]=1; df2[df2<=0]=0
# add nHostgene and nTFBS
df2=cbind(df2, nTFBS_gt5=ifelse(df$f06.TFBS>5,1,0))
#df2=cbind(df2, nHostgene_gt10=ifelse(df$f22.lenHostgeneMetaintron>0 & 1000 * df$f20.nHostgene / df$f22.lenHostgeneMetaintron >= 0.095, 1, 0))
df2[df2==0]=2; df2[df2==1]=0; df2[df2==2]=1  # 0<-->1  in order to get s10=(a+b+c)/(a+b+c+d)
#plot(hclust(dist.binary(df2,method=10), method='single'), labels=F)
DF2=dist.binary(df2,method=10)

## standarize the variable
df1=scale(df1, center = T)
df2=scale(df2, center = T) / sqrt(2)

df=cbind(df1,df2)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
distfunc2 <- function(x) as.dist((1 - cor(x))/2)
distfunc3 <- function(x) dist.binary(x,method=10)

cl.row <- hclustfunc(distfunc(df))
cl.col <- hclustfunc(distfunc2(df))
#plot(cl.row)
#plot(cl.col)
# extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
gr.row <- cutree(cl.row, 8)
gr.col <- cutree(cl.col, 8)

# require(RColorBrewer)
col1 <- brewer.pal(8, "Set1")
col2 <- brewer.pal(8, "Pastel1")

collist<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
my_palette <- colorRampPalette(collist)(n = 1000)

xmin=-5; xmax=5
df[df<xmin]=xmin; df[df>xmax]=xmax;

heatmap.2(as.matrix(df),
          lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1, 5, 1 ), lwid =c(1, 5 ),
          #cellnote = mat_data,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(9,4),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          cexCol=1,
          symbreak=F,
          scale='none',
          cexRow = 0.7,
          labRow = NA, #df0$freq, #[rownames(mat_data)],
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          #hclustfun=hclustfunc, distfun=distfunc,  
          Colv = as.dendrogram(cl.col),
          Rowv = as.dendrogram(cl.row),
          keysize = 0.8,
          key.title = NULL,
          key.xlab = NULL,
          key.ylab = NULL)

# df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
# df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
# df[df>0]=1;
# df=cbind(df, nEnhancer=apply(df, 1, sum))
# df=cbind(df, nHostgene=features$f20.nHostgene)
# df=cbind(df, cHostgene=cut(features$f20.nHostgene, breaks=c(-1,0,9,100000), labels=F))
# df=cbind(df, RPM=log10(features$f03.RPM), RPKM=log10(features$f02.RPKM))
# df=df[with(df, order(-nHostgene, -nEnhancer, -f09.chromHMM_brain, -f12.DNaseROADMAP, -f08.CAGEenhancer, -f07.P300, -f06.TFBS, -f15.HCNE)),]
# df=df[with(df, order(-f09.chromHMM_brain, -f08.CAGEenhancer, -nEnhancer, -f12.DNaseROADMAP, -f07.P300, -f06.TFBS, -f15.HCNE, nHostgene)),]
# 
# ## update version for progress report
# pdf("HiTNE.clustering.v3.pdf", width=2.5, height=8, paper='us')
# nheat=8
# par(mar=c(2,0,2,0), oma=c(3,3,3,3))
# layout(matrix(seq(nheat*2),nrow=2,ncol=nheat),widths=rep(.5,8),heights=c(3,0.3))
# 
# x=df$nHostgene
# xmin=0; xmax=200;
# x[x<xmin]=xmin; x[x>xmax]=xmax;
# collist<-c("#ffffff", "#000000")
# ColorRamp<-colorRampPalette(collist)(10000)
# ColorLevels<-seq(from=xmin, to=xmax, length=10000)
# ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(3,1,cex=2, las=2, labels="# of HiTNEs in host gene", hadj=1, padj = .5,tck=0, mgp=c(3,-.1,0))
# #axis(4,seq(0,length(x),5000),cex=2, labels=seq(0,length(x),5000)/1000, padj = .5,tck=-0.04, mgp=c(3,0.2,0))
# #mtext(paste0("HiTNE index (n=", nrow(df),")"), 4, line=2, cex=.8)
# axis(2,seq(0,length(x),5000),cex=2, labels=seq(0,length(x),5000)/1000, padj = .5,tck=-0.04, mgp=c(3,0.3,0))
# mtext(paste0("HiTNE index (n=", nrow(df),")"), 2, line=1, cex=.8)
# 
# barplot(-x, horiz=T, col=colorRampPalette(c('blue','red'))(100), lend=2, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# 
# image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,seq(0,10000,length=5),labels=seq(xmin,xmax,length=5),mgp=c(mgp=c(3,0.2,0)))
# 
# x=df$nEnhancer
# collist<-c("white", colorRampPalette(c("#D1E5F0","#053061"))(6))
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="nEnhancer",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(3,1,cex=2, las=2, labels="# of enhancer evidences", hadj=1, padj = .5,tck=0, mgp=c(3,-.1,0))
# #axis(4,seq(0,length(x),5000),cex=2, labels=seq(0,length(x),5000)/1000, padj = .5,tck=-0.04, mgp=c(3,0.2,0))
# #mtext(paste0("HiTNE index (n=", nrow(df),")"), 4, line=2, cex=.8)
# image(1:length(collist),1,matrix(data=1:length(collist),nrow=length(collist),ncol=1),col=collist, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,1:length(collist),labels=0:6,mgp=c(mgp=c(3,0.2,0)))
# 
# featureslist=c("f09.chromHMM_brain", "f08.CAGEenhancer", "f12.DNaseROADMAP", "f07.P300","f06.TFBS", "f15.HCNE");
# labelslist=c("chromHMM","CAGE","DNase","P300","TFBS","Conservation");
# colorslist=c("#ED6120", "#E44892", "#3370C0", "#644591", "#FF3851", "#2B3F87");
# for(i in 1:length(featureslist)){
#     x=df[,featureslist[i]]
#     collist<-c("white", colorslist[i])
#     image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
#     axis(3,1,cex=2, las=2, labels=labelslist[i], hadj=1, padj = .5,tck=0, mgp=c(3,-.1,0))
#     image(as.matrix(1),col=colorslist[i], xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# }
# 
# dev.off()
# 
# df=subset(df, nEnhancer>0)
# pdf("HiTNE.clustering.v4.pdf", width=2.5, height=8, paper='us')
# nheat=7
# par(mar=c(2,0,2,0), oma=c(3,3,3,3))
# layout(matrix(seq(nheat*2),nrow=2,ncol=nheat),widths=rep(.5,7),heights=c(3,0.3))
# 
# x=df$nEnhancer
# collist<-c(colorRampPalette(c("#D1E5F0","#053061"))(6))
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="nEnhancer",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(3,1,cex=2, las=2, labels="# of enhancer evidences", hadj=1, padj = .5,tck=0, mgp=c(3,-.1,0))
# #axis(4,seq(0,length(x),5000),cex=2, labels=seq(0,length(x),5000)/1000, padj = .5,tck=-0.04, mgp=c(3,0.2,0))
# #mtext(paste0("HiTNE index (n=", nrow(df),")"), 4, line=2, cex=.8)
# axis(2,seq(0,length(x),5000),cex=2, labels=seq(0,length(x),5000)/1000, padj = .5,tck=-0.04, mgp=c(3,0.3,0))
# mtext(paste0("HiTNE index (n=", nrow(df),")"), 2, line=1, cex=.8)
# image(1:length(collist),1,matrix(data=1:length(collist),nrow=length(collist),ncol=1),col=collist, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,1:length(collist),labels=1:6,mgp=c(mgp=c(3,0.2,0)))
# 
# featureslist=c("f09.chromHMM_brain", "f08.CAGEenhancer", "f12.DNaseROADMAP", "f07.P300","f06.TFBS", "f15.HCNE");
# labelslist=c("chromHMM","CAGE","DNase","P300","TFBS","Conservation");
# colorslist=c("#ED6120", "#E44892", "#3370C0", "#644591", "#FF3851", "#2B3F87");
# for(i in 1:length(featureslist)){
#     x=df[,featureslist[i]]
#     collist<-c("white", colorslist[i])
#     image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
#     axis(3,1,cex=2, las=2, labels=labelslist[i], hadj=1, padj = .5,tck=0, mgp=c(3,-.1,0))
#     image(as.matrix(1),col=colorslist[i], xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# }
# 
# dev.off()
# 
# 
# 
# df=df[with(df, order(-cHostgene,nEnhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f08.CAGEenhancer,f06.TFBS, f07.P300,f15.HCNE, RPM)),]
# pdf("HiTNE.clustering.pdf", width=8, height=8, paper='us')
# nheat=5
# par(mar=c(2,0,2,0), oma=c(3,3,3,3))
# layout(matrix(seq(nheat*2),nrow=2,ncol=nheat),widths=c(0.5,0.5,2,.5,.5),heights=c(3,0.3))
# 
# x=df$cHostgene
# collist<-c("#228B22","#FEBAAD","#FF4500")
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(1,1,cex=2, labels="# of HiTNE\nin host gene", padj = .5)
# mtext(paste0("HiTNEs (n=", nrow(df),")"), 2, line=1)
# image(1:length(collist),1,matrix(data=1:length(collist),nrow=length(collist),ncol=1),col=collist, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,1:length(collist),labels=c("0","[1,9]", ">=10"),mgp=c(mgp=c(3,0.2,0)))
# 
# x=df$nEnhancer
# collist<-c("white", colorRampPalette(c("#D1E5F0","#053061"))(6))
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), col=collist, useRaster=F, xlab="nEnhancer",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(1,1,cex=2, labels="# of enhancer\nevidences", padj = .5)
# image(1:length(collist),1,matrix(data=1:length(collist),nrow=length(collist),ncol=1),col=collist, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,1:length(collist),labels=0:6,mgp=c(mgp=c(3,0.2,0)))
# 
# x=df[,grep("^f", colnames(df))]
# x=x[,c("f09.chromHMM_brain", "f12.DNaseROADMAP", "f08.CAGEenhancer", "f06.TFBS", "f07.P300","f15.HCNE")]
# collist<-c(0,1)
# image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1:ncol(x),labels=c("chromHMM","DNase","CAGE","TFBS","P300","Conservation"), cex=.6)
# #text(1:ncol(x),y=-1,srt = -45, cex=.8, pos=1, labels=c("chromHMM","DNase","CAGE","TFBS","P300","Conservation"), xpd = TRUE)
# image(1:length(collist),1,matrix(data=1:length(collist),nrow=length(collist),ncol=1),col=collist, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,1:length(collist),labels=c("absent","present"),mgp=c(mgp=c(3,0.2,0)))
# 
# x=df$RPM
# xmin=-1.1; xmax=.5;
# x[x<xmin]=xmin; x[x>xmax]=xmax;
# collist<-c("#dddddd","#000000")
# ColorRamp<-colorRampPalette(collist)(10000)
# ColorLevels<-seq(from=xmin, to=xmax, length=10000)
# ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
# axis(1,1, labels=c("log10(RPM)"))
# image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,seq(0,10000,length=5),labels=seq(xmin,xmax,length=5),mgp=c(mgp=c(3,0.2,0)))
# 
# x=df$RPKM
# xmin=1.5; xmax=3;
# x[x<xmin]=xmin; x[x>xmax]=xmax;
# collist<-c("#dddddd","#000000")
# ColorRamp<-colorRampPalette(collist)(10000)
# ColorLevels<-seq(from=xmin, to=xmax, length=10000)
# ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1, labels=c("log10(RPKM)"))
# image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
# axis(1,seq(0,10000,length=4),labels=seq(xmin,xmax,length=4),mgp=c(mgp=c(3,0.2,0)))
# 
# 
# dev.off()
# 
# library(plyr); df=count(df, vars = colnames(df))
# df=df[with(df, order(-cHostgene,nEnhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f08.CAGEenhancer,f06.TFBS, f07.P300,f15.HCNE)),]
# 
# pdf("HiTNE.merged.clustering.pdf", height=0, width=1.38,paper='us')
# nheat=4
# par(mar=c(9,0,0,0), oma=c(2,3,2,2))
# layout(matrix(seq(nheat),nrow=1,ncol=nheat),widths=c(.5,.5,3,.5))
# 
# x=df$cHostgene
# collist<-c("#228B22","#FEBAAD","#FF4500")
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1,las=2,cex.axis = .9, labels="# of HiTNE in host gene", padj = .5)
# mtext(paste0("HiTNE clusters (n=", nrow(df),")"), 2, line=1)
# 
# x=df$nEnhancer
# collist<-c("white", colorRampPalette(c("#D1E5F0","#053061"))(6))
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=collist, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1,las=2,cex.axis = .9, labels="# of enhancer evidences", padj = .5)
# 
# x=df[,grep("^f[0-9]", colnames(df))]
# x=x[,c("f09.chromHMM_brain", "f12.DNaseROADMAP", "f08.CAGEenhancer", "f06.TFBS", "f07.P300","f15.HCNE")]
# collist<-c(0,1)
# image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=collist, asp=1, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1:ncol(x),las=2,cex.axis = .9, labels=c("chromHMM","DNase","CAGE","TFBS","P300","Conservation"))
# 
# x=df[,'freq']
# xmin=10; xmax=1000;
# x[x<xmin]=xmin; x[x>xmax]=xmax;
# collist<-c("#dddddd","#000000")
# ColorRamp<-colorRampPalette(collist)(10000)
# ColorLevels<-seq(from=xmin, to=xmax, length=10000)
# ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
# image(1,1:length(x), matrix(x, nrow=1, ncol=length(x)), asp=1, col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
# axis(1,1,las=2,cex.axis = .9, labels=c("# of HiTNEs"))
# 
# dev.off()


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