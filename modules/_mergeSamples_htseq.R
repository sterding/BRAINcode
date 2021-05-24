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
# 2021/1/5: remove the version .xx in the ENSID (as somehow different GENCODE version was used before and after batch 5 in BRAINcode. This might cause some genes lost in the follwing intersect() step
#fpkm$tracking_id = gsub("\\.\\d+$", "", fpkm$tracking_id)

for(i in 2:(n-1)){
    message(paste("[Merging file", args[i], "...]", i,"/",n-1, "Done"));
    df=read.table(args[i], header=F); colnames(df)=c("tracking_id", "FPKM");
    #df$tracking_id = gsub("\\.\\d+$", "", df$tracking_id)
    common=intersect(fpkm$tracking_id, df$tracking_id)
    fpkm=cbind(fpkm[match(common,fpkm$tracking_id),], FPKM=df[match(common,df$tracking_id), 'FPKM'])
}

# reassign column name
colnames(fpkm)[grep('FPKM',colnames(fpkm))]=gsub(".*run_output/([^/]*)/.*", "\\1", args[1:(n-1)])

# save data
write.table(fpkm, outputfile, sep="\t", quote = F, col.names = T, row.names = F)



