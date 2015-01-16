###########################################
# Rscript to check the consistency of RNAseq replicates
# Usage: Rscript $PATH/_replicate.R genes.fpkm.allSamples.uniq.xls genes.fpkm.allSamples.uniq.rep.pdf
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Nov-12
###########################################

args<-commandArgs(TRUE)

FPKMfile=args[1]  # either filename or stdin
outputfile=args[2]

# FPKMfile="genes.fpkm.HCILB.uniq.xls"; outputfile="genes.fpkm.HCILB.uniq.QC.pdf" 

message("loading data...")

fpkm=read.table(file(FPKMfile), header=T, check.names =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; fpkm=fpkm[,grep("FPKM",colnames(fpkm))]; colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))

replicated=gsub("(.*_.*_.*)_.*_.*","\\1", colnames(fpkm))[duplicated(gsub("(.*_.*_.*)_.*_.*","\\1", colnames(fpkm)))]
message(paste(length(replicated),"samples have found with replicates!"))

pdf(outputfile, width=6, height=6)
for(i in replicated){
    message(paste("plot for sample", i, "..."))
    ii=grep(i, colnames(fpkm))
    plot(fpkm[,ii] + 1e-5, log='xy',pch='.', cex=0.6, main=paste("Two replicates for",i), xlab=colnames(fpkm)[ii[1]], ylab=colnames(fpkm)[ii[2]])
    legend("topleft", paste("Pearson's r =", round(cor(fpkm[,ii])[2], 3)), bty='n')
}
dev.off()
