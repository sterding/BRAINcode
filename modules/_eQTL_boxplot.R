###############################################
## Rscript for generating boxplot of expression vs. genotype for gene-SNP pairs (e.g. from the eQTL output)
## Author: Xianjun Dong
## Date: 2016-May-26
## Version: 0.0
## Usage: Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_boxplot.R gene.snp.list path
###############################################
require(MatrixEQTL)

args<-commandArgs(TRUE)
GSfile=args[1]  # gene SNP table
path=ifelse(is.na(args[2]),getwd(),args[2]) 
setwd(path); 

# setwd("~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/"); load("data.RData"); load("genes.RData"); G="ENSG00000069275.12"; S="rs823116:205720483:G:A_G:A"; genesnp = read.table("final.cis.eQTL.d1e6.p1e-2.xls", header=T, stringsAsFactors =F)
# setwd("~/eRNAseq/HCILB_SNDA");load("data.RData"); load("genes.RData"); G="chr17_44218414_44219042"; S="rs17649553:43994648:C:T_C:T"; genesnp = read.table("final.cis.eQTL.d1e6.p1e-2.xls", header=T, stringsAsFactors =F)

message("# load data...")
######################
load("data.RData") # snps etc.
genesnp = read.table("final.cis.eQTL.d1e6.p1e-2.xls", header=T, stringsAsFactors =F)
if(file.exists("genes.RData")) load("genes.RData") else{
  residuals = read.table("expression.postSVA.xls", check.names=F)
  genes = SlicedData$new();
  genes$CreateFromMatrix(as.matrix(residuals))
  save(genes, file="genes.RData")
}

GS=read.table(GSfile, header = F, stringsAsFactors = F)
apply(GS, 1, function(gs) {
  G=gs[1]; S=gs[2];
  message(paste0("# making eQTL plot for ",S," and ",G," ..."))
  ######################
  RS=sub("([^:]*):.*","\\1", S) ## the part before the first ":" in the SNP ID
  REF = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][1]  ## get the REF allele
  ALT = strsplit(sub(".*_(.*)","\\1", S),":")[[1]][2]  ## get the ALT allele
  
  df=data.frame(expression=as.numeric(genes$FindRow(G)$row), SNP=as.numeric(snps$FindRow(S)$row), row.names = colnames(snps$FindRow(S)$row))
  ## write data to txt file
  write.table(df, file=paste("eQTLboxplot",G,gsub(":","..",S),"xls",sep="."), col.names = NA,quote=F, sep="\t")

  genesnp0=subset(genesnp, gene==G & SNP==S)
  p=signif(genesnp0$p.value, 3);
  
  df$SNP1=factor(df$SNP, levels=2:0)  ## in the current All.Matrix.txt, the number is number of REF allele (since we use REF in the --recode-allele) -- Ganqiang
  if(is.na(p) || length(p)==0) { # in case SNP:eQTL pair is not in the final.cis.eQTL.d1e6.p1e-2.xls file (e.g. those not or less significant pairs)
      # test=aov(expression~SNP1, df)
      # p=signif(summary(test)[[1]][["Pr(>F)"]][[1]],3)  
      # change to use lm(): 12/11/2018
      test = lm(expression~SNP1, df)
      p=signif(summary(test)[[1]][["Pr(>t)"]][[1]],3)  
  }
  pdf(file=paste("eQTLboxplot",G,gsub(":","..",S),"pdf",sep="."), width=6, height=4)
  par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
  bp=boxplot(expression~SNP1, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
  stripchart(expression~SNP1, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=1, add = TRUE) 
  title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
  mtext(paste0(c(REF,REF,ALT),c(REF,ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:3)
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
    mtext(c(paste0(c(REF,REF),c(REF,ALT),collapse = "/"), paste0(ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
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
    mtext(c(paste0(REF,REF), paste0(c(REF,ALT),c(ALT,ALT),collapse = "/")),  cex=0.7,side=1,line=1.5,at=1:2)
    mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  }
  
  mtext(paste("cis-eQTL for",G,"and",S), outer = TRUE, cex = 1.2)     
  
  ## without jitter and outliers (optional)
  ##############################

  par(mfrow=c(1,3), mar=c(4,4,2,1), oma = c(0, 0, 2, 0))
  bp=boxplot(expression~SNP1, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outline = F, outpch=NA)
  title(main=paste0("additive effect (p = ", p,")"), cex.main=1, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), cex=0.7,side=1,line=.5,at=1:3)
  mtext(paste0(c(REF,REF,ALT),c(REF,ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:3)
  mtext(paste0("N = ", bp$n),  cex=0.7,side=1,line=2.5,at=1:3)
  
  # with/without REF
  if(length(unique(df$SNP2))>1 & min(table(df$SNP2))>1) { #grouping factor must have exactly 2 levels  AND at least 2 values in each group
      p0=signif(t.test(expression~SNP2, df)$p.value,3)
      bp=boxplot(expression~SNP2, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
      title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
      mtext(c("with REF","without REF"), cex=0.7,side=1,line=0.5,at=1:2)
      mtext(c(paste0(c(REF,REF),c(REF,ALT),collapse = "/"), paste0(ALT,ALT)),  cex=0.7,side=1,line=1.5,at=1:2)
      mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  }
  
  # with/without ALT
  df$SNP3=ifelse(as.numeric(as.character(df$SNP))==2,0,1)  ## 0: without ALT; 1:with ALT
  if(length(unique(df$SNP3))>1 & min(table(df$SNP3))>1) {  #grouping factor must have exactly 2 levels AND at least 2 values in each group
      p0=signif(t.test(expression~SNP3, df)$p.value,3)
      bp=boxplot(expression~SNP3, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outline = F, outpch=NA)
      title(main=paste0("dominant effect (p = ",p0,")"), cex.main=1, line=0.5)
      mtext(c("without ALT","with ALT"), cex=0.7,side=1,line=0.5,at=1:2)
      mtext(c(paste0(REF,REF), paste0(c(REF,ALT),c(ALT,ALT),collapse = "/")),  cex=0.7,side=1,line=1.5,at=1:2)
      mtext(paste0("N = ", bp$n), cex=0.7,side=1,line=2.5,at=1:2)
  }
  
  mtext(paste("cis-eQTL for",G,"and",S), outer = TRUE, cex = 1.2)     
  
  dev.off() 
  
})