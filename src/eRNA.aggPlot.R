## script to analysis enhancer regions

setwd("~/eRNAseq")

mark="eRNA"
marks=c("RNAseq/RNAseq.BRAINCODE.HCILB_SNDA","CAGE/CAGE.FANTOM5.total.rev","CAGE/CAGE.FANTOM5.total.fwd","DNase/DNase.Roadmap.fBrain.pval","Histone/H3K27ac.Roadmap.SN.pval","Histone/H3K4me1.Roadmap.SN.pval", "TFBS/TFBS.ENCODE.all.count","Conservation/Conservation.phyloP46way")

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
output_dir="~/eRNAseq/";

## read the binned data
X=data.frame()
for(i in 1:length(marks))
{
  x0=read.table(paste0(output_dir, "externalData/", marks[i], ".bigwig.",mark,".1kbp.summit.or.mid.1kbp.100bins"), nrows =71469, check.names =F, stringsAsFactors=F, header=F); rownames(x0)=x0[,1]; x0=x0[,-1];
  if(ncol(X)>0) {X=cbind(X,x0);} else { X=x0;}
}

## read the classification of eRNA candidate
features = read.table("~/eRNAseq/eRNA.characterize.xls", nrows =71469, check.names =F, stringsAsFactors=F, header=T)
df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
df$class=ifelse(apply(df,1,sum)==0,1,ifelse(df$f12.DNaseROADMAP==0,ifelse(apply(df,1,sum)==1,2,ifelse(apply(df,1,sum)==2,3,4)),5))

# enhancers = subset(features, f06.TFBS>=5 | f07.P300>0 | f08.CAGEenhancer>0 | f09.chromHMM_brain>0 | f12.DNaseROADMAP>0 | f15.HCNE>0, select = c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
# enhancers$f06.TFBS=ifelse(enhancers$f06.TFBS>5,1,0)
# enhancers[enhancers>0]=1
# coreEnhancers = enhancers[apply(enhancers,1,sum)>=3,]
# 
# hclustfunc <- function(x) hclust(x, method="single")
# distfunc <- function(x) dist(x,method="euclidean")
# distfunc2 <- function(x) as.dist((1 - cor(x, method='spearman'))/2)
# distfunc3 <- function(x) dist.binary(x,method=10)
# x0=X[rownames(coreEnhancers),]
# cl <- hclustfunc(distfunc((x0)))
# ORDER1 <- cl$order
# 
# # randomly take 3398 rows from the one without any enhancer evidence (0)
# nonEnhancers = subset(features, f06.TFBS==0 & f07.P300==0 & f08.CAGEenhancer==0 & f09.chromHMM_brain==0 & f12.DNaseROADMAP==0 & f15.HCNE== 0)
# dim(nonEnhancers)
# nonEnhancers=nonEnhancers[sample.int(nrow(nonEnhancers),nrow(coreEnhancers)),]
# x0=X[rownames(nonEnhancers),]
# dim(x0)
# cl <- hclustfunc(distfunc((x0)))
# ORDER2 <- cl$order

MARKS=c("RNAseq","CAGE(-)","CAGE(+)","DNase","H3K4me1","H3K27ac","nTFBS","phyloP")
cols=c("red", 'magenta', 'darkmagenta','darkgreen','darkorange4','blue','purple','darkblue')
LOG =c(1, 1, 1,    1, 1,   1,    0, 0)
MIN=c(-2, 0,   0,   1,   -1, -0.5,   0, 0)
MAX=c( 0, 0.9, 0.6, 3.5,  1.5, 1,   5, 1)
ONLYMAX=c(0,0,0,0,0,0,0,0)

png("aggregation.eRNAs.class.all.png",width=1000, height=800)
par(mfrow=c(6,8), mar=c(.5,0,0,0), oma=c(2,1,1,1))
layout(matrix(seq(6*length(marks)),nrow=6,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=rep(1,6))
ymin=c(0.03,0.68,0.68,0,0.2,0.2,0,0.05)
ymax=c(0.1,1.45,1.45,8.5,2.2,2.3,10,0.3)
for(k in c(0:5))
{
  if(k==0) index=rownames(X) else { index=rownames(subset(df,class==k))}
  for(i in 1:length(marks))
  {
    x=X[index,(100*(i-1)+1):(i*100)]
    plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
    V=ncol(x)/2+.5
    abline(v=V, lty=2, lwd=.5)
  }
}
dev.off()

png("aggregation.eRNAs.class.four.png",width=1200, height=200)
pdf("aggregation.eRNAs.class.four.pdf",width=18, height=2, paper = 'usr')
par(mfrow=c(1,8), mar=c(.5,0,2,0), oma=c(2,1,1,1))
layout(matrix(seq(1*length(marks)),nrow=1,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=1)
ymin=c(0.03,0.68,0.68,0,0.2,0.2,0,0.05)
ymax=c(0.1,1.05,1.05,8.5,2.2,2.3,2,0.18)
for(i in 1:length(marks))
{
  # all
  x=X[,(100*(i-1)+1):(i*100)]
  plot(apply(x, 2, mean, trim =0.01), main=MARKS[i], type='l',col='black', ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
  
  # none
  x=X[rownames(subset(df, class==1)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='blue')
  
  # non-DNase
  x=X[rownames(subset(df, class>1 & class<5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='purple')

  # DNase
  x=X[rownames(subset(df, class==5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='red')
  
  V=ncol(x)/2+.5
  abline(v=V, lty=2, lwd=.5)
}
dev.off()
  
for(k in c(0:5))
{
  if(k==0) index=rownames(X) else { index=rownames(subset(df,class==k))}
  
  D1=apply(X[index,], 1, function(x) (sum(x[1:50])-sum(x[51:100]))/(sum(x[1:50])+sum(x[51:100])))  # RNAseq
  ORDER1=order(D1)
  
  message(paste("processing the class",k))
  
  png(paste0("aggregation.eRNAs.class",k,".all.png"),width=1000, height=800)
  par(mfcol=c(3,8), mar=c(.5,0,0,0), oma=c(2,1,1,1))
  layout(matrix(seq(3*length(marks)),nrow=3,ncol=length(marks)),widths=rep(0.5,length(marks)),heights=c(0.5,3,0.2))
  

  for(i in 1:length(marks))
  {
    x=X[index,(100*(i-1)+1):(i*100)]
    plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
    V=ncol(x)/2+.5
    abline(v=V, lty=2, lwd=.5)
    
    x=x[ORDER1, ];
    if(ONLYMAX[i]) x=apply(x,1,max)
    if(LOG[i]) x=log10(x+0.01)
    range(x)
    xmin=range(x)[1];xmax=range(x)[2]
    xmin=MIN[i];xmax=MAX[i];
    x[x<xmin]=xmin;  x[x>xmax]=xmax;
    collist<-c("white",cols[i])
    ColorRamp<-colorRampPalette(collist)(1000)
    ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*1000/(xmax-xmin)) : round( (max(x)-xmin)*1000/(xmax-xmin) )]
    ColorLevels<-seq(from=xmin, to=xmax, length=100)
    if(ONLYMAX[i]==0) image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
    if(ONLYMAX[i]==1) image(1,1:length(x),matrix(data=x,nrow=1, ncol=length(x)),col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
    abline(v=V, lty=2, lwd=.5)
    
    image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=colorRampPalette(collist)(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
    if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0), cex=0.5, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex=0.5, hadj=1)}
    text(50,1,gsub(".*/(.*)","\\1",MARKS[i]))
  }
  dev.off()
}


#pdf("aggregation.enhancers.SN.pdf", width=18, height=16)
png("aggregation.enhancers.SN.png",res=300, width=672, height=607)
par(mfcol=c(1,9), mar=c(.5,0,0,0), oma=c(2,1,1,1))
layout(matrix(seq(4*length(marks)),nrow=4,ncol=length(marks)),widths=rep(.5,length(marks)),heights=c(0.5,3,3,0.2))
#layout(matrix(seq(length(marks)),nrow=1,ncol=length(marks)),widths=rep(.5,length(marks)))


MARKS=c("RNAseq","CAGE.fwd","CAGE.rev","DNase","H3K4me1","H3K27ac","nTFBS","Conservation")
cols=c("red", 'darkmagenta','magenta', 'darkgreen','orange','blue','purple','darkblue')
LOG =c(1, 0, 0,    1, 1,   1,    0, 0)
MIN=c(-2, 0, 0, -1.2,-1.2, -1.5, 0, 0)
MAX=c(1,  2, 2,  1.2, 1.2,  1.5, 5, 1)
ONLYMAX=c(0,0,0,0,0,0,0,0)

for(i in 1:length(marks))
{
  x0=X[,(100*(i-1)+1):(i*100)]
  plot(apply(x0, 2, mean, trim =0.01), type='l',col=cols[i], xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  
  x=x0[rownames(coreEnhancers),]
  x=x[ORDER1, ];
  if(ONLYMAX[i]) x=apply(x,1,max)
  if(LOG[i]) x=log10(x+0.01)
  range(x)
  xmin=range(x)[1];xmax=range(x)[2]
  xmin=MIN[i];xmax=MAX[i];
  x[x<xmin]=xmin;  x[x>xmax]=xmax;
  #barplot(apply(x, 2,mean));
  # boxplot(log(x))
  collist<-c("white",cols[i])
  ColorRamp<-colorRampPalette(collist)(1000)
  ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*1000/(xmax-xmin)) : round( (max(x)-xmin)*1000/(xmax-xmin) )]
  ColorLevels<-seq(from=xmin, to=xmax, length=100)
  if(ONLYMAX[i]==0) image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  if(ONLYMAX[i]==1) image(1,1:length(x),matrix(data=x,nrow=1, ncol=length(x)),col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)

  # random non-enhancers
  x=x0[rownames(nonEnhancers),]
  x=x[ORDER2, ];
  if(ONLYMAX[i]) x=apply(x,1,max)
  if(LOG[i]) x=log10(x+0.01)
  range(x)
  xmin=range(x)[1];xmax=range(x)[2]
  xmin=MIN[i];xmax=MAX[i];
  x[x<xmin]=xmin;  x[x>xmax]=xmax;
  #barplot(apply(x, 2,mean));
  # boxplot(log(x))
  collist<-c("white",cols[i])
  ColorRamp<-colorRampPalette(collist)(1000)
  ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*1000/(xmax-xmin)) : round( (max(x)-xmin)*1000/(xmax-xmin) )]
  ColorLevels<-seq(from=xmin, to=xmax, length=100)
  if(ONLYMAX[i]==0) image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  if(ONLYMAX[i]==1) image(1,1:length(x),matrix(data=x,nrow=1, ncol=length(x)),col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  
  image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=colorRampPalette(collist)(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0), cex=0.5, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex=0.5, hadj=1)}
  text(50,1,gsub(".*/(.*)","\\1",MARKS[i]))
}

dev.off()