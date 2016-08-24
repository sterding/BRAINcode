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
S=args[1]  # SNP
G=args[2]  # Gene

#G="chr10_61632545_61633312"; S="rs1684902:61633127:A:G_A:G"

message("# load data...")
load("data.RData")

RS=sub("([^:]*):.*","\\1", S) ## the part before the first ":" in the SNP ID
REF = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][1]  ## get the REF allele
ALT = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][2]  ## get the ALT allele

message("# making eQTL plot ...")
######################
genesnp = read.table("final.cis.eQTL.d1e6.p1e-2.xls", header=T, stringsAsFactors =F)
if(file.exists("genes.RData")) load("genes.RData") else{
  residuals = read.table("expression.postSVA.xls", check.names=F)
  genes = SlicedData$new();
  genes$CreateFromMatrix(as.matrix(residuals))
  save(genes, file="genes.RData")
}

# one file for all genes (only the most significant SNP per gene)
pdf(paste("eQTLplot",G,S,"pdf",sep="."), width=6, height=4)

par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
genesnp0=subset(genesnp, gene==G & SNP==S)
p=signif(genesnp0$p.value, 3);

df=data.frame(expression=as.numeric(genes$FindRow(G)$row), SNP=as.numeric(snps$FindRow(S)$row))
df$SNP1=factor(df$SNP, levels=2:0)  ## in the current All.Matrix.txt, the number is number of REF allele (since we use REF in the --recode-allele) -- Ganqiang
if(is.na(p)) p=signif(t.test(expression~SNP1, df)$p.value,3)  # in case SNP:eQTL pair is not in the final.cis.eQTL.d1e6.p1e-2.xls file (e.g. those less significant)
bp=boxplot(expression~SNP1, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
stripchart(expression~SNP1, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE) 
title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
mtext(paste0("(",c(REF,REF,ALT),c(REF,ALT,ALT),")"),  cex=0.7,side=1,line=1.5,at=1:3)
mtext(paste0("N = ", bp$n),  cex=0.7,side=1,line=2.5,at=1:3)

# with/without REF
df$SNP2=ifelse(as.numeric(as.character(df$SNP))==0,0,1)  ## 0: without REF; 1:with REF
if(length(unique(df$SNP2))>1 & min(table(df$SNP2))>1) { #grouping factor must have exactly 2 levels  AND at least 2 values in each group
    df$SNP2=factor(df$SNP2, levels=1:0)
    p0=signif(t.test(expression~SNP2, df)$p.value,3)
    bp=boxplot(expression~SNP2, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
    stripchart(expression~SNP2, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE)
    title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
    mtext(c("with REF","without REF"), cex=0.7,side=1,line=0.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=1.5,at=1:2)
}

# with/without ALT
df$SNP3=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
if(length(unique(df$SNP3))>1 & min(table(df$SNP3))>1) {  #grouping factor must have exactly 2 levels AND at least 2 values in each group
    df$SNP3=factor(df$SNP3, levels=0:1)
    p0=signif(t.test(expression~SNP3, df)$p.value,3)
    bp=boxplot(expression~SNP3, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
    stripchart(expression~SNP3, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE)
    title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
    mtext(c("without ALT","with ALT"), cex=0.7,side=1,line=0.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=1.5,at=1:2)
}

mtext(paste("cis-eQTL for",G,"and",RS), outer = TRUE, cex = 1.2)     

## without jitter and outliers (optional)
##############################

par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
bp=boxplot(expression~SNP1, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outline = F, outpch=NA)
title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
mtext(paste0("(",c(REF,REF,ALT),c(REF,ALT,ALT),")"),  cex=0.7,side=1,line=1.5,at=1:3)
mtext(paste0("N = ", bp$n),  cex=0.7,side=1,line=2.5,at=1:3)

# with/without REF
if(length(unique(df$SNP2))>1 & min(table(df$SNP2))>1) { #grouping factor must have exactly 2 levels  AND at least 2 values in each group
    p0=signif(t.test(expression~SNP2, df)$p.value,3)
    bp=boxplot(expression~SNP2, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
    title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
    mtext(c("with REF","without REF"), cex=0.7,side=1,line=0.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=1.5,at=1:2)
}

# with/without ALT
df$SNP3=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
if(length(unique(df$SNP3))>1 & min(table(df$SNP3))>1) {  #grouping factor must have exactly 2 levels AND at least 2 values in each group
    p0=signif(t.test(expression~SNP3, df)$p.value,3)
    bp=boxplot(expression~SNP3, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
    title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
    mtext(c("without ALT","with ALT"), cex=0.7,side=1,line=0.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=1.5,at=1:2)
}

mtext(paste("cis-eQTL for",G,"and",RS), outer = TRUE, cex = 1.2)     

dev.off() 
# 
# ## sample with heterzegous SNPs
# write.table(t(snps$FindRow(S)$row),file = paste0("SNP.",RS), quote=F, sep="\t", col.names = NA)
