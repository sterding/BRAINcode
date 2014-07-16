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

fpkm=read.table(args[1], header=T);
fpkm=subset(fpkm, select=c('tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length'))

for(i in 1:(n-1)){
    message(paste("[Merging file", args[i], "...] %", round(100*i/(n-1), 1), "Done"));
    df=read.table(args[i], header=T)
    fpkm=cbind(fpkm, FPKM=df[match(fpkm$tracking_id, df$tracking_id), 'FPKM'])
}

# reassign column name
colnames(fpkm)[grep('FPKM',colnames(fpkm))]=paste("FPKM", gsub(".*run_output/([^/]*)/.*", "\\1", args[1:(n-1)]), sep=".")

# save data
write.table(fpkm, outputfile, sep="\t", quote = F, col.names = T, row.names = F)



