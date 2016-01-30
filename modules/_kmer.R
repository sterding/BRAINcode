###########################################
# Rscript to generate the kmer distance plot
# Usage: Rscript $PATH/_kmer.R mergedHCILB_k9.matrix.txt mergedHCILB_k9.matrix.pdf
# Reference: https://www.genevestigator.com/userdocs/manual/qc.html
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Aug-8
###########################################

args<-commandArgs(TRUE)

matrix=args[1]  # 'mergedHCILB_k9.matrix.txt'
outputfile=args[2]
if(is.na(outputfile)) outputfile=paste(matrix,"pdf",sep=".")

message("reading the matrix data...")

N=read.table(matrix, header=F,nrows=1)$V1
nms=read.table(matrix, header=F, skip=1, nrows=N)$V1
df=data.matrix(read.table(matrix,  fill = TRUE, skip=N+1, col.names = 1:(N-1)))
upper=rbind(cbind(0,t(df)),0); upper[is.na(upper)]=0
lower=rbind(0,cbind(df,0)); lower[is.na(lower)]=0
df=lower+upper; diag(df)=NA; dim(df)
kmer=apply(df,1,function(x) median(x,na.rm=T))
names(kmer)=nms

pdf(outputfile, width=7, height=7)
hist(kmer, breaks=100, xlab="Median k-mer distance", ylab="Number of samples")
legend("topright", paste(names(kmer[kmer>0.001]), kmer[kmer>0.001], sep=": "), bty='n', cex=0.8)
dev.off()
