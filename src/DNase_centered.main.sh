cd ~/eRNAseq/externalData/DNase

## call narrow peaks from merged uniform signal of individual samples [SELECTED as final solution]
grep DNase macs2signal.list | cut -f2 | awk '{print "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"$1}' | xargs -n 1 -P 8 wget -b
mkdir download; mv E* download/
bsub -q big -n 1 -M 10000 bigWigMerge download/E*DNase.pval.signal.bigwig merged.DNase.pval.signal.bg
bsub -q big -n 1 -M 10000 bedGraphToBigWig merged.DNase.pval.signal.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size merged.DNase.pval.signal.bigwig
module load macs2/2.1
bsub -q big -n 1 -M 10000 macs2 bdgpeakcall -i merged.DNase.pval.signal.bg --min-length 100 -o merged.DNase.pval.signal.p0.01.halfsamples.peaks --cutoff 53  # at least half samples with p=0.01 (2*53/2)
ln -fs merged.DNase.pval.signal.p0.01.peaks merged.DNase.pval.signal.peaks
intersectBed -a merged.DNase.pval.signal.peaks -b ../../eRNA.bed -u | wc -l
#6129

# ## merge narrow peaks called from individual samples [DISCARDED eventually due to two ]
# grep DNase macs2signal.list | cut -f2 | sed 's/pval.signal.bigwig/macs2.narrowPeak.gz/g' | awk '{print "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"$1}' | xargs -n 1 -P 8 wget -b
# zcat *DNase.macs2.narrowPeak.gz | sortBed > merged.DNase.macs2.narrowPeak
# awk '{OFS="\t"; $10=$2+$10;print}' merged.DNase.macs2.narrowPeak | mergeBed -i - -d -50 -c 4,5,6,7,8,9,10 -o count,max,distinct,max,max,max,median| awk '{OFS="\t"; $10=$10-$2; print}' > merged.DNase.macs2.narrowPeak.merged
# intersectBed -a merged.DNase.macs2.narrowPeak.merged -b ../../eRNA.bed -u | awk '$4>=40' | wc -l 
# # 3075

# alternative: use the 638304 enhancer-being DNase peaks defined in 10 Roadmap "brain" samples (http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation)
awk '($3-$2)>100' regions_enh_merged.brain.bed | while read chr start end rest
do
  l=$(expr $end - $start)
  bigWigSummary merged.DNase.pval.signal.bigwig $chr $start $end $l | awk -vchr=$chr -vstart=$start -vend=$end '{OFS="\t"; for(i=1;i<=NF;i++) if($i!="n/a" && $i>max) {imax=i;max=$i}}END{print chr, start, end, ".", 0, ".", max, -1, -1, imax-1}'
done > regions_enh_merged.brain.narrowPeak

wc -l regions_enh_merged.brain.narrowPeak

## script to generate data to draw aggregation plot for HITNE

#============================================================
# get bin signal
#============================================================
cd ~/eRNAseq

intersectBed -a externalData/DNase/merged.DNase.pval.signal.peaks -b eRNA.bed -u | awk '{OFS="\t"; mid=$2+$10; print $1,mid-1000,mid+1000, $4}' > eRNA.with.mergedDNase.peak.summit.1kbp.bed

for i in externalData/*/*.bigwig;
do
    toBinRegionsOnBigwig.sh $i eRNA.with.mergedDNase.peak.summit.1kbp.bed 100 > $i.eRNA.with.mergedDNase.peak.summit.1kbp.100bins &
done

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.aggPlot.R


## script to analysis enhancer regions

setwd("~/eRNAseq")

mark="eRNA"
marks=c("RNAseq/RNAseq.BRAINCODE.HCILB_SNDA","CAGE/CAGE.FANTOM5.total.rev","CAGE/CAGE.FANTOM5.total.fwd","DNase/merged.DNase.pval.signal","Histone/H3K27ac.Roadmap.SN.pval","Histone/H3K4me1.Roadmap.SN.pval", "TFBS/TFBS.ENCODE.all.count","Conservation/Conservation.phyloP46way")

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
output_dir="~/eRNAseq/";

X=data.frame()
for(i in 1:length(marks))
{
  x0=read.table(paste0(output_dir, "externalData/", marks[i], ".bigwig.",mark,".with.mergedDNase.peak.summit.1kbp.100bins"), check.names =F, stringsAsFactors=F, header=F); 
  rownames(x0)=x0[,1]; x0=x0[,-1];
  if(ncol(X)>0) {X=cbind(X,x0);} else { X=x0;}
}

hclustfunc <- function(x) hclust(x, method="single")
distfunc <- function(x) dist(x,method="euclidean")
distfunc2 <- function(x) as.dist((1 - cor(x, method='spearman'))/2)
distfunc3 <- function(x) dist.binary(x,method=10)

cl <- hclustfunc(distfunc((X[,1:100])))
#cl <- hclustfunc(distfunc2(t(X)))
ORDER1 <- cl$order

D1=apply(X, 1, function(x) (sum(x[1:50])-sum(x[51:100]))/(sum(x[1:50])+sum(x[51:100])))  # RNAseq
D2=apply(X, 1, function(x) (sum(x[101:150])-sum(x[151:200]))/(sum(x[101:150])+sum(x[151:200])))  # CAGE rev
D3=apply(X, 1, function(x) (sum(x[201:250])-sum(x[251:300]))/(sum(x[201:250])+sum(x[251:300])))  # CAGE fwd
ORDER2=order(D1)

pdf("aggregation.enhancers.all.1kbp.orderbyD.pdf", width=8, height=6, paper='us')
#png("aggregation.enhancers.all.1kbp.orderbyD.png", width=800, height=600)
#png("aggregation.enhancers.all.1kbp.orderbyCOR.png", width=800, height=600)
par(mfcol=c(1,9), mar=c(.5,0,0,0), oma=c(2,1,1,1))
layout(matrix(seq(3*length(marks)),nrow=3,ncol=length(marks)),widths=rep(.5,length(marks)),heights=c(0.4,3,0.2))

MARKS=c("RNAseq","CAGE(-)","CAGE(+)","DNase","H3K4me1","H3K27ac","nTFBS","phyloP")
cols=c("red", 'magenta', 'darkmagenta','darkgreen','darkorange4','blue','purple','darkblue')
LOG =c(1, 1, 1,    1, 1,   1,    0, 0)
MIN=c(-2, 0,   0,   1,   -1, -0.5,   0, 0)
MAX=c( 0, 0.9, 0.6, 3.5,  1.5, 1,   5, 1)
ONLYMAX=c(0,0,0,0,0,0,0,0)

for(i in 1:length(marks))
{
  x=X[,(100*(i-1)+1):(i*100)]
  plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
  V=ncol(x)/2+.5
  abline(v=V, lty=2, lwd=.5)
  
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
  abline(v=V, lty=2, lwd=.5)
  
  image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=colorRampPalette(collist)(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0),  cex.axis=1, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex.axis=1,  hadj=1)}
  text(50,1,MARKS[i], cex=1.5)
}
dev.off()


# png("aggregation.enhancers.all.500bp.png", width=800, height=600)
# par(mfcol=c(1,9), mar=c(.5,0,0,0), oma=c(2,1,1,1))
# layout(matrix(seq(3*length(marks)),nrow=3,ncol=length(marks)),widths=rep(.5,length(marks)),heights=c(0.5,3,0.2))
# 
# MARKS=c("RNAseq","CAGE.rev","CAGE.fwd","DNase","H3K4me1","H3K27ac","nTFBS","Conservation")
# cols=c("red", 'magenta', 'darkmagenta','darkgreen','orange','blue','purple','darkblue')
# LOG =c(1, 0, 0,    1, 1,   1,    0, 0)
# MIN=c(-2, 1, 1, 1,  -1, -1, 0, 0)
# MAX=c(1,3.5,3.5, 3.5,  3,  3,  5, 1)
# ONLYMAX=c(0,0,0,0,0,0,0,0)
# 
# for(i in 1:length(marks))
# {
#   x=X[,(100*(i-1)+1):(i*100)]
#   x=x[,26:75]
#   plot(apply(x, 2, mean, trim =0.01), type='l',col=cols[i], xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
#   V=ncol(x)/2+.5
#   abline(v=V, lty=2, lwd=.5)
#   
#   x=x[ORDER1, ];
#   if(ONLYMAX[i]) x=apply(x,1,max)
#   if(LOG[i]) x=log10(x+0.01)
#   range(x)
#   xmin=range(x)[1];xmax=range(x)[2]
#   xmin=MIN[i];xmax=MAX[i];
#   x[x<xmin]=xmin;  x[x>xmax]=xmax;
#   #barplot(apply(x, 2,mean));
#   # boxplot(log(x))
#   collist<-c("white",cols[i])
#   ColorRamp<-colorRampPalette(collist)(1000)
#   ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*1000/(xmax-xmin)) : round( (max(x)-xmin)*1000/(xmax-xmin) )]
#   ColorLevels<-seq(from=xmin, to=xmax, length=100)
#   if(ONLYMAX[i]==0) image(1:ncol(x),1:nrow(x), as.matrix(t(x)), col=ColorRamp_ex, useRaster=F, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
#   if(ONLYMAX[i]==1) image(1,1:length(x),matrix(data=x,nrow=1, ncol=length(x)),col=ColorRamp_ex, useRaster=T, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n", frame.plot=T)
#   abline(v=V, lty=2, lwd=.5)
#   
#   image(1:length(ColorLevels),1,matrix(data=1:length(ColorLevels),nrow=length(ColorLevels),ncol=1),col=colorRampPalette(collist)(100), xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
#   if(i==1) {axis(1, 1,labels=" low", mgp=c(3,0.2,0),  cex.axis=1.5, hadj = 0);axis(1, 100,labels="high ", mgp=c(3,0.2,0), cex.axis=1.5,  hadj=1)}
#   text(50,1,gsub(".*/(.*)","\\1",MARKS[i]), cex=2)
# }
# dev.off()


## how many eRNA candidates overlap with a DNase peak
intersectBed -b externalData/DNase/merged.DNase.pval.signal.peaks -a eRNA.bed -u > eRNA.wDNase.bed
#5669

## how many of them has a bimodel around the DNase summit
cat eRNA.with.mergedDNase.peak.summit.1kbp.bed | while read chr start end name rest
do
  bigWigSummary externalData/RNAseq/RNAseq.BRAINCODE.HCILB_SNDA.bigwig $chr $start $end 100 | _filter.awk -vW=5 -vtype=mean | _wave_detector.awk -vchr=$chr -vstart=$start -vend=$end | awk -vstart=$start -vend=$end '{OFS="\t"; mid=(start+end)/2; min=($7<$9)?$7:$9; if(mid<$3 && mid>$2 && (min/($8+0.001))>=2) print}'
done > eRNA.with.mergedDNase.peak.with.svs.bed

wc -l eRNA.with.mergedDNase.peak.with.svs.bed
# 2039

awk '{OFS="\t"; min=($7<$9)?$7:$9; max=($7<$9)?$9:$7; if((min/($8+0.001))>=5 && max<(min*2)) print}' eRNA.with.mergedDNase.peak.with.svs.bed | wc -l
# 350
awk '{OFS="\t"; min=($7<$9)?$7:$9; max=($7<$9)?$9:$7; if((min/($8+0.001))>=5 && max<(min*2)) print}' eRNA.with.mergedDNase.peak.with.svs.bed | intersectBed -a stdin -b externalData/CAGE/permissive_enhancers.bed -u | wc -l
# 33

## Note: it seems not many HiTNEs are bimodel around the peak of DNase