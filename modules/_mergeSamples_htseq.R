###########################################
# Rscript to merge samples'd expression data from cufflinks
# Usage: Rscript $PIPELINE_PATH/_mergeSamples.R `ls ../../run_output/*/uniq/genes.fpkm_tracking` mergedFPKM.tab
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-May-22
###########################################

args<-commandArgs(TRUE)

n=length(args)
outputfile=args[n]

fpkm=read.table(args[1], header=F);
colnames(fpkm)=c("tracking_id", "FPKM")

for(i in 2:(n-1)){
    message(paste("[Merging file", args[i], "...] %", round(100*i/(n-1), 1), "Done"));
    df=read.table(args[i], header=F); colnames(df)=c("tracking_id", "FPKM");
    fpkm=cbind(fpkm, FPKM=df[match(fpkm$tracking_id, df$tracking_id), 'FPKM'])
}

# reassign column name
colnames(fpkm)[grep('FPKM',colnames(fpkm))]=paste("FPKM", gsub(".*run_output/([^/]*)/.*", "\\1", args[1:(n-1)]), sep=".")

# save data
write.table(fpkm, outputfile, sep="\t", quote = F, col.names = T, row.names = F)



