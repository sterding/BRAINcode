###########################################
# Rscript to infer target gene for eRNA by expression correlation
# Usage: Rscript $0 eRNA.expression.tab genes.fpkm.allSamples.uniq.xls output_filename
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Sep-12
###########################################

args<-commandArgs(TRUE)

df1=read.table(args[1], header=T)
df2=read.table(args[2], header=T)  
output_filename=args[3]
method=args[4] # "pearson"

# df1=read.table("eRNAfinal.allsamples.RPKM.tab",header=T); df2=read.table("genes.fpkm.allSamples.uniq.xls",header=T)

COR_CUTOFF = 0.7

rownames(df1)=df1[,1]; df1=df1[,-1]; rownames(df2)=df2[,1]; df2=df2[,-1];
colnames(df1) = gsub("FPKM.","",colnames(df1),ignore.case = T); colnames(df2) = gsub("FPKM.","",colnames(df2),ignore.case = T)
df1=df1[,order(colnames(df1))]; df2=df2[,order(colnames(df2))]

# filter out eRNAs or genes with low mean values
df1=df1[rowMeans(df1)>=0.005,]; df2=df2[rowMeans(df2)>=0.005,]
# to check if the distribution is skewed
# hist(log(rowMeans(df2)), breaks=100)

# transpose
df1=t(df1); df2=t(df2)

## use pearson on log transformed (because pearson can be misleading on Highly skewed variables, see http://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data)
#ddcor1=cor(log(0.005+df1),log(0.005+df2), method="pearson")
# or, use spearman
ddcor=cor(df1,df2, method="spearman")

# filter eRNA-gene pairs based on COR_CUTOFF
ddcor=ddcor[apply(ddcor,1,max) > COR_CUTOFF,]
ddcor.filted = do.call(rbind,apply(cbind(1:nrow(ddcor),ddcor),1,function(x) {y=x[-1];if(sum(y>COR_CUTOFF)>0) {cbind(rownames(ddcor)[x[1]],colnames(ddcor)[y>COR_CUTOFF],round(y[y>COR_CUTOFF],3))}}))

# write to file (for visualzation in Circos)
write.table(ddcor.filted, paste(output_filename,".rho.tab",sep="."), col.names = c("eRNA","gene","rho"), row.names=F,sep="\t", quote=F)

# filter eRNA-gene pairs based on Pearson's r > COR_CUTOFF
ddcor=cor(log(0.005+df1),log(0.005+df2), method="pearson")
ddcor=ddcor[apply(ddcor,1,max) > COR_CUTOFF,]
ddcor.filted = do.call(rbind,apply(cbind(1:nrow(ddcor),ddcor),1,function(x) {y=x[-1];if(sum(y>COR_CUTOFF)>0) {cbind(rownames(ddcor)[x[1]],colnames(ddcor)[y>COR_CUTOFF],round(y[y>COR_CUTOFF],3))}}))

# write to file (for visualzation in Circos)
write.table(ddcor.filted, paste(output_filename,".pcc.tab",sep="."), col.names = c("eRNA","gene","pcc"), row.names=F,sep="\t", quote=F)


