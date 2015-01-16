###########################################
# Rscript to plot the expression vs. genotyping for significant eQTL SNP-gene pairs
# Usage: Rscript $PIPELINE_PATH/_eQTL.plot.R output.cis.txt snp.txt expression.txt cov.txt
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-May-26
###########################################

args<-commandArgs(TRUE)

output_file_name=args[1]  # output of Matrix-eQTL 
SNP_file_name=args[2]
expression_file_name=args[3]
covariates_file_name=args[4]

pCutoff=as.numeric(args[5])
if(is.na(pCutoff)) pCutoff=0.05

print(pCutoff)

genesnp = read.table(output_file_name, header=T);
genesnp = subset(genesnp, FDR<=pCutoff);

snps = read.table(SNP_file_name, header=T); snps=t(data.frame(snps[,-1], row.names=snps[,1]));  
ge   = read.table(expression_file_name, header=T); ge=t(data.frame(ge[,-1], row.names=ge[,1])); 
cov  = read.table(covariates_file_name, header=T); cov=t(data.frame(cov[,-1], row.names=cov[,1]));

pdf(paste(output_file_name, "expression.vs.genotype","pdf", sep="."), width=10, height=4)
par(mfrow=c(1,3), mar=c(4,4,2,2), oma = c(0, 0, 2, 0))
for(i in 1:nrow(genesnp)){
    s=as.character(genesnp$SNP[i]); g=as.character(genesnp$gene[i]);
    print(paste(i, s, g))
    df=as.data.frame(cbind(cov, expression=ge[,g], SNP=snps[,s]))
    df$SNP=factor(df$SNP, levels=0:2)
    #fit=lm(expression ~ condition + age + sex + SNP, data=df)
    #df=cbind(res=residuals(fit), df)
    df$condition=ifelse(df$condition==0,"HC", ifelse(df$condition==1, "ILB", "PD"))
    df$condition=factor(df$condition, levels=c("HC","ILB","PD"))
    
    #bp=boxplot((1+expression)~condition*SNP, data=df, ylab="Expression (FPKM)", col=c('#2ca25f','#feb24c','#ff0000'),log="y", at=c(1,2,3,5,6,7,9,10,11), xaxt='n', main="")
    bp=boxplot(expression~condition*SNP, data=df, ylab="Rank Normalized Gene Expression", col=c('#2ca25f','#feb24c','#ff0000'),at=c(1,2,3,5,6,7,9,10,11), xaxt='n', main="")
    title(main=paste(paste("SNP:",s),paste("GENE:",g),sep="\n"), cex.main=0.8, line=0.5)
    axis(side=1,at=c(2,6,10),cex = 1, labels=c("0/0","0/1","1/1"), cex=2)
    # count of persons in each category
    mtext(c("N:", bp$n), side=1,line=2,at=c(0, 1,2,3,5,6,7,9,10,11),cex = .8)

    bp=boxplot(expression~SNP, data=df, ylab="Rank Normalized Gene Expression", xaxt='n', main="")
    title(main=paste(paste("SNP:",s),paste("GENE:",g),sep="\n"), cex.main=0.8, line=0.5)
    axis(side=1,cex = 1, at=1:3, labels=c("0/0","0/1","1/1"), cex=2)
    mtext(c("N:", bp$n), side=1,line=2,at=c(0:3),cex = .8)

    df$SNP=ifelse(as.numeric(as.character(df$SNP))==0,0,1)
    df$SNP=factor(df$SNP, levels=0:1)
    bp=boxplot(expression~condition*SNP, data=df, ylab="Rank Normalized Gene Expression", col=c('#2ca25f','#feb24c','#ff0000'),at=c(1,2,3,5,6,7), xaxt='n', main="")
    title(main=paste(paste("SNP:",s),paste("GENE:",g),sep="\n"), cex.main=0.8, line=0.5)
    axis(side=1,at=c(2,6),cex = 1, labels=c("w/ allele","w/o allele"), cex=2)
    mtext(c("N:", bp$n), side=1,line=2,at=c(0,1,2,3,5,6,7),cex = .8)
}
dev.off()