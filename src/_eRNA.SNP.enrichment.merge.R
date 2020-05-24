###########################################
# Rscript to merge samples'd expression data from cufflinks
# Usage: Rscript $PIPELINE_PATH/_eRNA.SNP.enrichment.merge.R eRNA.SNP.enrichments.SNAP HCILB_SNDA HC_nonNeuron HC_PY
# Note: can choose between only the enriched ones (eRNA.SNP.enrichments.SNAP.xls) vs. all (eRNA.SNP.full.SNAP.xls)
# Author: Xianjun Dong
# Version: 0.1
# Date: 2016-12-21
###########################################
setwd("~/eRNAseq")

args<-commandArgs(TRUE)

n=length(args)
whichfile=args[1]

fpkm=read.table(paste0(args[2],"/", whichfile, ".xls"), header=T, stringsAsFactors=F);
fpkm$type=ifelse(fpkm$type=="HTNE", paste0(fpkm$type,".",args[2]), fpkm$type)

for(i in 3:n){
    message(paste("[Merging file", args[i], "...] %", round(100*(i-1)/n, 1), "Done"));
    df=read.table(paste0(args[i],"/", whichfile, ".xls"), header=T, stringsAsFactors=F)
    df=subset(df, type=="HTNE")
    df$type=paste0(df$type,".",args[i])
    
    fpkm=rbind(fpkm, df)
}

# save data
#options(digit = 3);
write.table(format(fpkm, digit=3), paste0(whichfile, ".merge.xls"), sep="\t", na="", row.names = F, quote =F)

# HTNE only and in wide format
# =================================
## elegant version
# library(dplyr)
# library(reshape2)
# library(data.table)
# fpkm=read.table("eRNA.SNP.enrichment.merge.xls", sep="\t", header=T, stringsAsFactors=F);
# fpkm = fpkm %>% filter(grepl("HTNE",type)) %>% as.data.table() %>%
#   data.table::dcast(Disease_or_Trait ~ type, value.var=c('observed','all','pvalue','OR'), sep=".") # dcast in data.table can do multiple varibles (Ref: http://stackoverflow.com/a/29068012/951718)
# head(fpkm)

fpkm=read.table(paste0(args[2],"/", whichfile, ".xls"), header=T, stringsAsFactors=F);
fpkm=subset(fpkm, type=="HTNE", select=-c(type))
names(fpkm)=ifelse(names(fpkm)=="Disease_or_Trait",names(fpkm), paste0(names(fpkm),".",args[2]))

for(i in 3:n){
  message(paste("[Merging file", args[i], "...] %", round(100*(i-1)/n, 1), "Done"));
  df=read.table(paste0(args[i],"/", whichfile, ".xls"), header=T, stringsAsFactors=F)
  df=subset(df, type=="HTNE", select=-c(type))
  #fpkm = merge(fpkm, df, by="Disease_or_Trait", all=T, suffixes =c("1",paste0(".",args[i])))  # suffixes only works for the columns with the same names
  fpkm = merge(fpkm, df, by="Disease_or_Trait", all=T)  
  names(fpkm)[(ncol(fpkm)-ncol(df)+2):ncol(fpkm)] = paste0(names(fpkm)[(ncol(fpkm)-ncol(df)+2):ncol(fpkm)],".",args[i])
}

# save data
#options(digit = 5, scientific=T);
write.table(fpkm, paste0(whichfile, ".merge.wide.xls"), sep="\t", row.names = F, quote =F)

# filter with at least one cell type with p<=0.01/1037
write.table(fpkm[apply(select(fpkm, contains('pvalue')),1,min, na.rm=T)<=0.01/1037,], paste0(whichfile, ".merge.wide.qvalue0.01.xls"), sep="\t", row.names = F, quote =F)

# filter with at least one cell type with p<=0.05/1037
write.table(fpkm[apply(select(fpkm, contains('pvalue')),1,min, na.rm=T)<=0.05/1037,], paste0(whichfile, ".merge.wide.qvalue0.05.xls"), sep="\t", row.names = F, quote =F)