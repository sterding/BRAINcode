## Rscript to calculate the consistence of HTNE expression across the samples
## Usage: Rscript  /data/neurogen/pipeline/RNAseq/src/_HTNE.consistency.R HC_nonNeuron
## Author: Xianjun Dong
## Date: Oc-23-2015

### 2: distribution of random background
args<-commandArgs(TRUE)

SAMPLE_GROUP=args[1]

# samplelist is pre-generated in _combine_bigwig.sh script
samplelist=read.table(paste0("~/neurogen/rnaseq_PD/results/merged/samplelist.",SAMPLE_GROUP), stringsAsFactors =F)$V1

pdf("background.RNAseq.cummulative.plot.pdf")
# plot the cummulative plot for merged track
df=read.table(paste("~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized", SAMPLE_GROUP, "bw.rdbg",sep="."), header=F)[,2] # mean RPM (mean0 from bigWigAverageOverBed)
Fn=ecdf(log10(df+0.0001))
plot(Fn, verticals = TRUE, do.points = FALSE, main=paste("trimmedmean.uniq.normalized", SAMPLE_GROUP,sep="."), xlab="average RPM", ylab="cummulative percentage (approx. 1-p)")
Fn=ecdf(log10(df[df>0]))
plot(Fn, verticals = TRUE, do.points = FALSE, add=T, col='blue')
legend("bottomright", c("all values (pseudocount: 0.0001)","non-zero values"), col=c('black','blue'), lty=1)
# define the cutoff with p=0.05 significance
inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
abline(h=0.95, v=g(0.95), col='red', lty=2, lwd=1)
points(g(0.95), 0.95, col='red', pch=19)
text(g(0.95), 0.95, round(g(0.95),2), cex=5, adj=c(0,1))
dev.off()

EXP=data.frame(); PV=data.frame(); QV=data.frame(); 
for(i in samplelist){
  print(i)
  # read background
  df=read.table(paste0("~/neurogen/rnaseq_PD/run_output/",i,"/uniq/accepted_hits.normalized2.bw.rdbg"), header=F)[,2] # mean RPM (mean0 from bigWigAverageOverBed)
  Fn=ecdf(df)
  
#   # plot the cummulative plot
#   plot(Fn, verticals = TRUE, do.points = FALSE, main=i, ylim=c(0.99, 1), xlab="average RPM", ylab="cummulative percentage (approx. 1-p)")
#   inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
#   abline(h=0.999, v=g(0.999), col='red', lty=2, lwd=1)
#   points(g(0.999), 0.999, col='red', pch=19)
#   text(g(0.999), 0.999, round(g(0.999),2), cex=5, adj=c(0,1))
  
  # read expression
  expression=read.table(paste0("~/neurogen/rnaseq_PD/run_output/",i,"/uniq/accepted_hits.normalized2.bw.eRNA.meanRPM"), header=F)
  pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));
  qvalue=as.numeric(format(p.adjust(pvalue, "BH"), digits=3));

  # merge
  if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; expression[,2]=qvalue; QV=expression; }
  else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); QV=cbind(QV, qvalue); }
}

colnames(EXP)=c("locus",samplelist); colnames(PV)=c("locus",samplelist); colnames(QV)=c("locus",samplelist);
write.table(EXP, "eRNA.meanRPM.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(PV,  "eRNA.pvalues.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(QV,  "eRNA.qvalues.xls", col.names=T, row.names=F, sep="\t", quote=F)

#rM=rowMeans(QV[,-1]<=0.05)
# write.table(EXP[rM>0.25,], paste("eRNA",SAMPLE_GROUP,"meanRPM.xls", sep="."), col.names=T, row.names=F, sep="\t", quote=F)
# write.table(PV[rM>0.25,],  paste("eRNA",SAMPLE_GROUP,"pvalue.xls", sep="."), col.names=T, row.names=F, sep="\t", quote=F)
# write.table(QV[rM>0.25,],  paste("eRNA",SAMPLE_GROUP,"qvalue.xls", sep="."), col.names=T, row.names=F, sep="\t", quote=F)
# 
# pdf(paste("eRNA",SAMPLE_GROUP,"qvalue.hist.pdf", sep="."), width=8, height=6)
# h=hist(rM, breaks=80, xlim=c(0,1), main="",xlab=paste("Percentage of",SAMPLE_GROUP, "samples (N=", length(samplelist),"with q-value <= 0.05"), ylab="Count of HTNEs", freq=T)
# abline(v=0.250, lty=2, col='red')
# legend('topright', c(bquote(.(sum(rM>0.25)) ~ "HTNEs"), expression("with q-value" <= "0.05"), "in at least 25% of samples"),  bty='n', text.col='red', cex=1.5)
# dev.off()