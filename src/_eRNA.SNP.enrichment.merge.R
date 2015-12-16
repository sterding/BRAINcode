###########################################
# Rscript to merge samples'd expression data from cufflinks
# Usage: Rscript $PIPELINE_PATH/_eRNA.SNP.enrichment.merge.R HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY 
# Author: Xianjun Dong
# Version: 0.0
# Date: 2015-11-15
###########################################

args<-commandArgs(TRUE)

n=length(args)

fpkm=read.table(paste0(args[1],"/eRNA.SNP.enrichments.SNAP.xls"), header=T, stringsAsFactors=F);
fpkm$type=ifelse(fpkm$type=="HTNE", paste0(fpkm$type,".",args[1]), fpkm$type)

for(i in 2:n){
    message(paste("[Merging file", args[i], "...] %", round(100*i/(n-1), 1), "Done"));
    df=read.table(paste0(args[i],"/eRNA.SNP.enrichments.SNAP.xls"), header=T, stringsAsFactors=F)
    df=subset(df, type=="HTNE")
    df$type=paste0(df$type,".",args[i])
    
    fpkm=rbind(fpkm, df)
}

# save data
#options(digit = 3);
write.table(format(fpkm, digit=3), "eRNA.SNP.enrichment.merge.xls", sep="\t", na="", row.names = F, quote =F)


# HTNE only and in wide format
# =================================
fpkm=read.table(paste0(args[1],"/eRNA.SNP.enrichments.SNAP.xls"), header=T, stringsAsFactors=F);
fpkm=subset(fpkm, type=="HTNE", select=-c(type))
names(fpkm)=ifelse(names(fpkm)=="Disease_or_Trait",names(fpkm), paste0(names(fpkm),".",args[1]))

for(i in 2:n){
  message(paste("[Merging file", args[i], "...] %", round(100*i/(n-1), 1), "Done"));
  df=read.table(paste0(args[i],"/eRNA.SNP.enrichments.SNAP.xls"), header=T, stringsAsFactors=F)
  df=subset(df, type=="HTNE", select=-c(type))
  fpkm = merge(fpkm, df, by="Disease_or_Trait", all=T, suffixes =c("",paste0(".",args[i])))
}

# save data
#options(digit = 5, scientific=T);
write.table(format(fpkm, digit=3), "eRNA.SNP.enrichment.merge.wide.xls", sep="\t", na="", row.names = F, quote =F)




