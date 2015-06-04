## script to analysis enhancer regions

mark="eRNA"
#============================================================
# draw aggregation plot (to run in R console)
#============================================================
output_dir="~/neurogen/rnaseq_PD/results/eRNA/externalData";
pdf("aggregation.enhancers.SN.pdf", width=18, height=16)
par(mfcol=c(1,9))

# HITNEs RPM
x=read.table(paste0(output_dir, "/RNAseq/eRNA.f03.RPM.txt"), header=F)
rownames(x)=x[,1]; colnames(x) = c("RPM"); x=x[,-1, drop=F];
x=log10(x);
hclust(dist(x))

features = read.table(paste0(output_dir, "/RNAseq/eRNA.characterize.xls"), header=T)
enhancers = subset(features, f06.TFBS>=5 | f07.P300>0 | f08.CAGEenhancer>0 | f09.chromHMM_brain>0 | f12.DNaseROADMAP>0 | f15.HCNE>0)
dim(enhancers)

ORDER = rownames(features)[order(features$f03.RPM, decreasing =T)]

marks=c("RNAseq/RNAseq","CAGE/CAGE.fwd","CAGE/CAGE.rev","DNase/E081-DNase.pval.signal","Histone/E074-H3K27ac.pval.signal","Histone/E074-H3K4me1.pval.signal","TFBS/TFBS","Conservation/Conservation")
MARKS=c("RNAseq","CAGE.fwd","CAGE.rev","DNase","H3K27ac","H3K4me1","TFBS","Conservation")
cols=c("red", 'magenta','magenta', 'darkgreen','orange','blue','purple','darkblue')
LOG =c(1, 0, 0,    1, 1,   1,    0, 0)
MIN=c(-2, 0, 0, -1.2,-1.2, -1.2, 0, 0)
MAX=c(1,  30, 30,  1.2, 1.2, 1.2,  5, 1)
ONLYMAX=c(0,1,1,0,0,0,0,0)

par(mar=c(1,0,0,0), oma=c(3,3,3,3))
layout(matrix(seq(2*length(marks)),nrow=2,ncol=length(marks)),widths=rep(.5,length(marks)),heights=c(3,0.2))
#layout(matrix(seq(length(marks)),nrow=1,ncol=length(marks)),widths=rep(.5,length(marks)))

for(i in 1:length(marks))
{
  x0=read.table(paste0(output_dir, "/", marks[i], ".bigwig.",mark,".1kbp.summit.or.mid.1kbp.100bins"), nrows =71469, check.names =F, stringsAsFactors=F, header=F); rownames(x0)=x0[,1]; x0=x0[,-1];
  #x=x[rownames(enhancers),]
  #ORDER=order(apply(x0,1,mean)); 
  x=x0[ORDER, ];
  if(ONLYMAX[i]) x=apply(x,1,sum)
  if(LOG[i]) x=log10(x+0.01)
  range(x)
  xmin=MIN[i];xmax=MAX[i];
  x[x<xmin]=xmin;  x[x>xmax]=xmax;
  #barplot(apply(x, 2,mean));
  # boxplot(log(x))
  collist<-c("white",cols[i])
  ColorRamp<-colorRampPalette(collist)(1000)
  ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*1000/(xmax-xmin)) : round( (max(x)-xmin)*1000/(xmax-xmin) )]
  ColorLevels<-seq(from=xmin, to=xmax, length=100)
  if(ONLYMAX[i]==0) image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  if(ONLYMAX[i]==1) image(1,1:length(x),matrix(data=x,nrow=1, ncol=length(x)),col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=colorRampPalette(collist)(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0), cex=0.5, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex=0.5, hadj=1)}
  text(50,1,gsub(".*/(.*)","\\1",MARKS[i]))
}