args<-commandArgs(TRUE)

df=read.table(args[1], header=F, stringsAsFactors = F); colnames(df) = c('Trait','count')
rownames(df)=df[,1]; df=df[,-1,drop=F]

for(i in args[2:length(args)]){
  message(i);
  d=try(read.table(i, header=F, stringsAsFactors = F)); 
  if(inherits(d, 'try-error')) {df=cbind(df, 0); next;}
  colnames(d) = c('Trait','count'); rownames(d)=d[,1]
  df=cbind(df, d[rownames(df),2])
  df[is.na(df)]=0
}
dim(df)
# empirical p-value = proportion of simulations in which the number of eQTLs exceeds the observed number
# ref: http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000888#s4
df= cbind(count=df[,1], pvalue=apply(df[,-1]>df[,1],1,mean))

write.table(df, file="dSNP.eSNP.enrichment.xls", quote = F, sep="\t")