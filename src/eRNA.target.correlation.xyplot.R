args<-commandArgs(TRUE)
inputfile=args[1]
outputfile=args[2]

df=read.table(inputfile, header=T)
rownames(df)=df[,1]; df=df[,-1]; df=t(df); colnames(df)=c("eRNA","target_gene"); df=as.data.frame(df)
df1=log(0.005+df)
reg1 <- lm(target_gene ~ eRNA, data=df1)
pdf(paste(outputfile,"pdf",sep='.'))
par(cex=.8)
plot(target_gene ~ eRNA, data=df1, main=outputfile)
abline(reg1)
legend("topleft",c(paste("Pearson's r:",round(cor(df1$target_gene,df1$eRNA, method="pearson"),3)), paste("Spearman's rho:",round(cor(df$target_gene,df$eRNA, method="spearman"),3))))
dev.off()
