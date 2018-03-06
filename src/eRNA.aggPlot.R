## script to analysis enhancer regions

setwd("~/eRNAseq/HCILB_SNDA")

#============================================================
# draw aggregation plot (to run in R console)
#============================================================

mark="eRNA"
marks=c("RNAseq/RNAseq.BRAINCODE.HCILB_SNDA",
        "DNase/merged.DNase.pval.signal",  # all Roadmap brain samples
        "Histone/H3k27ac.Roadmap.brainMerged", # all Roadmap brain samples
        "Histone/H3k4me1.Roadmap.brainMerged", # all Roadmap brain samples
        "Histone/H3k4me3.Roadmap.brainMerged", # all Roadmap brain samples
        "TFBS/TFBS.ENCODE.all.count",
        "CAGE/CAGE.FANTOM5.total.rev",
        "CAGE/CAGE.FANTOM5.total.fwd",
        "Conservation/Conservation.phyloP46way")

## read the binned data
X=data.frame()
for(i in 1:length(marks))
{
  binfile = paste0("~/neurogen/external_download/externalData/", marks[i], ".bigwig.",mark,".1kbp.summit.or.mid.1kbp.100bins")
  message(paste("Loading:",binfile))
  
  x0=read.table(binfile, check.names =F, stringsAsFactors=F, header=F); 
  rownames(x0)=x0[,1]; x0=x0[,-1];
  if(ncol(X)>0) {X=cbind(X,x0);} else { X=x0;}
}

## read the classification of eRNA candidate
features = read.table("eRNA.characterize.xls", check.names =F, stringsAsFactors=F, header=T)
library('dplyr')
df=features %>% select(f12.DNaseROADMAP, f09.chromHMM_brain, f07.P300, f06.TFBS, f08.CAGEenhancer, f15.HCNE) 
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
df$class=ifelse(apply(df,1,sum)==0,1,ifelse(df$f12.DNaseROADMAP==0,ifelse(apply(df,1,sum)==1,2,ifelse(apply(df,1,sum)==2,3,4)),5))
df$subclass=apply(df[,1:6],1,paste0,collapse="")

MARKS=c("RNAseq","DNase","H3K27ac","H3K4me1","H3K4me3","nTFBS","CAGE(-)","CAGE(+)","phyloP")
cols=c("red", 'darkgreen','darkorange4','darkorange4','blue','purple','magenta', 'darkmagenta','darkblue')

## to show aggregation of each class in separate rows
## -------------------------------------
# png("aggregation.eRNAs.class.all.png",width=1000, height=800)
# par(mfrow=c(6,8), mar=c(.5,0,0,0), oma=c(2,1,1,1))
# layout(matrix(seq(6*length(marks)),nrow=6,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=rep(1,6))
# ymin=c(0.03,0.68,0.68,0,0.2,0.2,0,0.05)
# ymax=c(0.1,1.45,1.45,8.5,2.2,2.3,10,0.3)
# for(k in c(0:5))
# {
#   if(k==0) index=rownames(X) else { index=rownames(subset(df,class==k))}
#   for(i in 1:length(marks))
#   {
#     x=X[index,(100*(i-1)+1):(i*100)]
#     plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
#     V=ncol(x)/2+.5
#     abline(v=V, lty=2, lwd=.5)
#   }
# }
# dev.off()

message("
        ## to show aggregation of each class in one row  (Fig. 3d)
        ## -------------------------------------
        ")
pdf("aggregation.eRNAs.class.four.pdf",width=18, height=2, paper = 'usr')
par(mfrow=c(1,9), mar=c(.5,0,2,0), oma=c(2,1,1,1), adj=c(0,0))
layout(matrix(seq(1*length(marks)),nrow=1,ncol=length(marks), byrow=1),widths=rep(1,length(marks)),heights=1)
ymin=c(0.03, 18,  1, 1,  1, 0.0, 0.8, 0.8, 0.03)
ymax=c(0.12, 180,10,13,13, 2.5, 1.4, 1.4, 0.2)
for(i in 1:length(marks))
{
  # all
  x=X[,(100*(i-1)+1):(i*100)]
  plot(apply(x, 2, mean, trim =0.01), main=MARKS[i], type='l',col='black', ylim=c(ymin[i],ymax[i]),xlab="",ylab="",xaxt="n",las=2, cex.axis=1,tck=0.04, frame.plot=T, mgp=c(3,-2,0))
  
  # class III (no mark)
  x=X[rownames(subset(df, class==1)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ffcccc')
  
  # class II (no DNase, but at least one other mark)
  x=X[rownames(subset(df, class>1 & class<5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ff9966')
  
  # class I (DNase)
  x=X[rownames(subset(df, class==5)),(100*(i-1)+1):(i*100)]
  points(x=1:100, y=apply(x, 2, mean, trim =0.01), type='l',col='#ff3333')
  
  V=ncol(x)/2+.5
  abline(v=V, lty=2, lwd=.5)
}

dev.off()


message("
        ## heatmaps of subclasses (Fig. 3e)
        ## -------------------------------------
        ")
for(subc in c('111110','111100','110001','011100','000000'))
{
  index=with(df,subclass==subc)
  XX=X[index,]
  
  if(subc=='000000'){
    XX=XX[order(apply(XX[,101:200],1,max)),]  # tune DNase signal a bit
    XX=XX[sample(dim(XX)[1]/2,1000),]
  }
  
  if(subc=='011100'){
    XX=XX[order(apply(XX[,101:200],1,max)),] # tune DNase signal a bit
    XX=XX[sample(dim(XX)[1]*0.9,min(1000,dim(XX)[1]*0.9)),]
  }
  
  D1=apply(XX, 1, function(x) (sum(x[1:50])-sum(x[51:100]))/(sum(x[1:50])+sum(x[51:100])))  # RNAseq
  ORDER1=order(D1)
  
  MARKS=c("RNAseq","DNase","H3K27ac","H3K4me1","H3K4me3","nTFBS","CAGE(-)","CAGE(+)","phyloP")
  cols=c("red", 'darkgreen','darkorange4','blue','lightblue','purple','magenta', 'darkmagenta','darkblue')
  
  LOG =c(0,    1,   0,   0,    0,   0, 0,   0,  0)
  MIN =c(0.04, 1.2,  2,    2,   2,   0, 1,   1,  0.1)
  MAX =c(0.1, 3.0,  15,   15,  10,  5,  3.5, 3.5, 0.8)
  if(subc=='110001') {MAX[2]=2.0; MIN[2]=1.1;}
  if(subc=='011100') {MAX[2]=4.0; MIN[2]=1.5;}
  
  message(paste("processing the subclass",subc))
  png(paste0("heatmap.eRNAs.subclass.",subc,".png"),width=1200, height=200)
  par(mfrow=c(1,length(marks)), mar=c(0,0,0,0), oma=c(0.1,0.1,0.1,0.1))
  layout(matrix(seq(1*length(marks)),nrow=1,ncol=length(marks), byrow=F),widths=rep(0.5,length(marks)),heights=3)
  for(i in 1:length(marks))
  {
    x=XX[,(100*(i-1)+1):(i*100)]
    
    #if(grepl("CAGE",marks[i])) x=x[,26:75]
    
    # histone of x
    #apply(x,2,range)
    
    x=x[ORDER1, ];
    if(LOG[i]) x=log10(x+0.01)
    
    # heatmap
    xmin=MIN[i];xmax=MAX[i];
    x[x<xmin]=xmin;  x[x>xmax]=xmax;
    collist<-c("white",cols[i])
    ColorRamp<-colorRampPalette(collist)(1000)
    image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
    V=ncol(x)/2+.5; abline(v=V, col=1,lty=2, lwd=1)
  }
  
  dev.off()
  
}

message("
        ## heatmaps of three classes (Fig. 3e)
        ## -------------------------------------
        ")
for(k in 1:4)
{
  if(k==1) index=rownames(subset(df,class==5))
  if(k==2) index=rownames(subset(df,class>1 & class<5))
  if(k==3) index=rownames(subset(df,class==1))
  if(k==4) index=rowSums(df[,1:6])>=3  # TNE with at least 3 suppored features
  
  XX=X[index,]
  
  D1=apply(XX, 1, function(x) (sum(x[1:50])-sum(x[51:100]))/(sum(x[1:50])+sum(x[51:100])))  # RNAseq
  ORDER1=order(D1)
  
  MARKS=c("RNAseq","DNase","H3K27ac","H3K4me1","H3K4me3","nTFBS","CAGE(-)","CAGE(+)","phyloP")
  cols=c("red", 'darkgreen','darkorange4','blue','lightblue','purple','magenta', 'darkmagenta','darkblue')
  
  LOG =c(0,    1,   0,   0,    0,   0, 0,   0,  0)
  MIN =c(0.04, 1.2,  2,    2,   2,   0, 1,   1,  0.1)
  MAX =c(0.1, 3.0,  15,   15,  10,  5,  3.5, 3.5, 0.8)
  
  message(paste("processing the class",k))
  png(paste0("heatmap.eRNAs.class.",k,".",paste(MAX,collapse="-"),".png"),width=800, height=600)
  
  par(mfrow=c(3,length(marks)), mar=c(.5,0,1,0), oma=c(2,1,1,1))
  layout(matrix(seq(3*length(marks)),nrow=3,ncol=length(marks), byrow=F),widths=rep(0.5,length(marks)),heights=c(.5,3,0.2))
  
  for(i in 1:length(marks))
  {
    x=XX[,(100*(i-1)+1):(i*100)]
    
    x=x[ORDER1, ];
    if(LOG[i]) x=log10(x+0.01)
    
    # aggregation
    plot(apply(x, 2, mean, trim =0.01), main="", type='l',col=cols[i], xlab="",ylab="",xaxt="n",las=2, cex.axis=1, tck=0.04, frame.plot=T, mgp=c(3,-2,0))
    V=ncol(x)/2+.5
    abline(v=V, lty=2, col=1, lwd=.5)
    
    # heatmap
    xmin=MIN[i];xmax=MAX[i];
    x[x<xmin]=xmin;  x[x>xmax]=xmax;
    collist<-c("white",cols[i])
    ColorRamp<-colorRampPalette(collist)(1000)
    image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n",frame.plot=T)
    V=ncol(x)/2+.5; abline(v=V, lty=2, col=1, lwd=.5)
    
    # color legend
    image(1:100,1,matrix(data=1:100,nrow=100,ncol=1),col=colorRampPalette(collist)(100), useRaster=F, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
    if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0), cex=0.5, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex=0.5, hadj=1)}
    text(50,1,MARKS[i], cex=1.5)
  }
  mtext(paste(paste(c("LOG:",LOG),collapse = ";"),paste(c("MIN:",MIN),collapse = ";"),paste(c("MAX:",MAX),collapse = ";")), outer = T, line=-.2, cex = .8)
  dev.off()
  
}