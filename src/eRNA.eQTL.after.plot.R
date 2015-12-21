###############################################
## Rscript for post-eQTL analysis (after running eRNA.eQTL.R)
## Author: Xianjun Dong
## Date: 2015-Dec-2
## Version: 0.0
## Usage: Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.plot.R rs1684902:61633127:A:G_A:G chr10_61632545_61633312
###############################################
#setwd("~/neurogen/rnaseq_PD/results/eQTL/eRNA/")
require(MatrixEQTL)

args<-commandArgs(TRUE)
S=args[1]  
G=args[2]

#G="chr10_61632545_61633312"; S="rs1684902:61633127:A:G_A:G"

message("# load data...")
load("data.RData")

RS=sub("([^:]*):.*","\\1", S) ## the part before the first ":" in the SNP ID
REF = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][1]  ## get the REF allele
ALT = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][2]  ## get the ALT allele

message("# making eQTL plot ...")
######################
genesnp = read.table("final.cis.eQTL.xls", header=T, stringsAsFactors =F)
residuals = read.table("expression.postSVA.xls")
genes = SlicedData$new();
genes$CreateFromMatrix(as.matrix(residuals))

# one file for all genes (only the most significant SNP per gene)
pdf(paste("eQTLplot",G,S,"mostsignificant.pdf",sep="."), width=8, height=8)
par(mfrow=c(1,2), mar=c(4,4,2,2), oma = c(0, 0, 2, 0))

  genesnp0=subset(genesnp, gene==G & SNP==S)
  p=signif(genesnp0$p.value, 3);
  
  df=data.frame(expression=as.numeric(genes$FindRow(G)$row), SNP=as.numeric(snps$FindRow(S)$row))
  df$SNP=factor(df$SNP, levels=2:0)  ## in the current All.Matrix.txt, the number is number of REF allele (since we use REF in the --recode-allele) -- Ganqiang
  bp=boxplot(expression~SNP, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
  stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE) 
  title(main=paste0("additive effect (pvalue=", p,")"), cex.main=0.8, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), side=1,line=0,at=1:3)
  mtext(paste0("(",c(REF,REF,ALT),c(REF,ALT,ALT),")"), side=1,line=1,at=1:3)
  mtext(paste0("N = ", bp$n), side=1,line=2,at=1:3)
  
  df$SNP=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
  df$SNP=factor(df$SNP, levels=0:1)
  p0=signif(t.test(expression~SNP, df)$p.value,3)
  bp=boxplot(expression~SNP, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
  stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE)
  title(main=paste0("dominant effect (pvalue=",p0,")"), cex.main=0.8, line=0.5)
  mtext(c("without ALT","with ALT"), side=1,line=0,at=1:2)
  mtext(paste0("N = ", bp$n), side=1,line=1,at=1:2)
  
  mtext(paste("cis-eQTL for",G,"and",RS), outer = TRUE, cex = 1.5)        
dev.off() 