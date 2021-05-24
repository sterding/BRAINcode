###########################################
# Rscript to merge samples'd expression data from cufflinks
# Usage: Rscript $PIPELINE_PATH/_mergeSamples.R `ls ../../run_output/*/uniq/genes.fpkm_tracking` mergedFPKM.tab
# Author: Xianjun Dong
# Version: 0.2
# Date: 2021-Jan-05
###########################################

args<-commandArgs(TRUE)

n=length(args)
outputfile=args[n]

fpkm=read.table(args[1], header=T, stringsAsFactors=T);
#fpkm=subset(fpkm, select=c('tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 'locus', 'length'))
fpkm=subset(fpkm, select=c('tracking_id', "FPKM")) # only fpkm 2018/05/23
# 2021/1/5: remove the version .xx in the ENSID (as somehow different GENCODE version was used before and after batch 5 in BRAINcode. This might cause some genes lost in the follwing intersect() step
#fpkm$tracking_id = gsub("\\.\\d+$", "", fpkm$tracking_id)

for(i in 2:(n-1)){
  message(paste("[Merging file", args[i], "...]", i,"/",n-1, "Done"));
  df=read.table(args[i], header=T, stringsAsFactors=T)
  #df$tracking_id = gsub("\\.\\d+$", "", df$tracking_id)
    ## only include the genes with FPKM_status==OK 
    #df=subset(df, FPKM_status=="OK")
  common=intersect(fpkm$tracking_id, df$tracking_id)
  fpkm=cbind(fpkm[match(common,fpkm$tracking_id),], FPKM=df[match(common,df$tracking_id), 'FPKM'])
}

# reassign column name
#colnames(fpkm)[grep('FPKM',colnames(fpkm))]=paste("FPKM", gsub(".*run_output/([^/]*)/.*", "\\1", args[1:(n-1)]), sep=".")
colnames(fpkm)[grep('FPKM',colnames(fpkm))]=gsub(".*run_output/([^/]*)/.*", "\\1", args[1:(n-1)])  # remove prefix FPKM 2018/03/05

# save data
write.table(fpkm, outputfile, sep="\t", quote = F, na="", row.names = F)



