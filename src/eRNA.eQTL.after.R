###############################################
## Rscript for post-eQTL analysis (after running eRNA.eQTL.R)
## Author: Xianjun Dong
## Date: 2015-Dec-2
## Version: 0.0
## Usage: bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.R
###############################################
#setwd("~/neurogen/rnaseq_PD/results/eQTL/eRNA/")
require(MatrixEQTL)
if(!require(qvalue)) {source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); require(qvalue);}

# load data
load("data.RData")

## 1. merge permutation result
observedP=read.table('final.cis.min.pv.gene.txt')
name=sort(rownames(observedP))
observedP=observedP[name,, drop=F]
min.pv.gene=read.table(paste0("permutations/permutation",1,".txt"))
for(i in list.files("permutations", pattern="permutation.*.txt", full.names=T)){
  message(i)
  min.pv.gene=cbind(min.pv.gene,read.table(i))
}
min.pv.gene=min.pv.gene[name, ]

#For each gene, empirical P value is set to the rank of observedP in the permutation list divided by X
# ref: https://www.biostars.org/p/13791/
empiricalP=apply(min.pv.gene <= observedP[,1], 1, mean)  # get the percenage of min(p) <= observedP among the 1000 permutations. 
#Note: the idea is similar as the empPvals() function in qvalue:  The p-values are calculated as the proportion of values from stat0 (or null p) that are greater (or smaller) than or equal to that from stat (or observed p).

#Use empirical P values as input to qvalue to estimate FDR for a given adjust P value cutoff.
q=qvalue(p = empiricalP)$qvalue

cat("# The estimate of the proportion of true alternative tests  (which is 3.5% in this case -- very low):", 1-qvalue(p = empiricalP)$pi0);

## now, we can get a list of all significantly associated SNP-gene pairs
#get the p-value correspinding to the eGene empirical p-value for eGenes at the q-value = 0.05 threshold
# range of min.pv.gene is not [0,1], which lead to "ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method" for following code
#cutoff = apply(min.pv.gene, 1, function(x) max(x[qvalue(x, lambda = seq(0,0.999,0.001))$qvalue<=0.05]))
cutoff = apply(min.pv.gene, 1, function(x) max(0, x[p.adjust(x, method='fdr')<=0.05]))

## or, to be simple, simply take the 5% percentile of emperical pvalue as cutoff
#cutoff = apply(min.pv.gene, 1, quantile, probs=0.05)

genesnp=read.table("final.cis.eQTL.xls", header=T, stringsAsFactors = FALSE)
genesnp=cbind(genesnp, SNPgene.cutoff=cutoff[match(genesnp$gene,names(cutoff))])
# filter on the eGene
genesnp=subset(genesnp, gene %in% names(q[q<=0.05]))
# filter on the SNP-gene cutoff
genesnp=subset(genesnp, p.value <= SNPgene.cutoff)

write.table(genesnp, file="final.cis.eQTL.FDR.05.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)

eGenes=aggregate(SNP~gene, data=data.frame(genesnp), FUN=length)
eGenes=cbind(eGenes, p.value=empiricalP[eGenes$gene], q.value=q[eGenes$gene])
write.table(eGenes[order(eGenes$gene, -eGenes$SNP),], "final.cis.eGene.qvalue.cutoff0.05.txt", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)


message("# making eQTL plot ...")
######################
genesnp = read.table("final.cis.eQTL.FDR.05.xls", header=T, stringsAsFactors =F)
residuals = read.table("expression.postSVA.xls")
genes = SlicedData$new();
genes$CreateFromMatrix(as.matrix(residuals))

# one file for all genes (only the most significant SNP per gene)
pdf("eQTLplot.mostsignificant.pdf", width=8, height=8)
par(mfrow=c(1,2), mar=c(4,4,2,2), oma = c(0, 0, 2, 0))
for(g in unique(genesnp$gene))
{
  genesnp0=subset(genesnp, gene==g)
  genesnp0 = genesnp0[with(genesnp0, order(FDR)), ]
  message(g);
  i=1
  s=as.character(genesnp0$SNP[i]);
  s2=sub("([a-zA-Z].*):(.*):(.*):(.*)","\\1", s) 
  g=genesnp0$gene[i];
  g_symbol=ifelse(is.null(genesnp0$gene_symbol[i]), genesnp0$gene[i], genesnp0$gene_symbol[i]);
  p=signif(genesnp0$p.value[i], 3);
  print(paste(i, s, s2, g_symbol, p))
  
  df=data.frame(expression=as.numeric(genes$FindRow(g)$row), SNP=as.numeric(snps$FindRow(s)$row))
  df$SNP=factor(df$SNP, levels=0:2)
  bp=boxplot(expression~SNP, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
  stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE) 
  title(main=paste0("additive effect (pvalue=", p,")"), cex.main=0.8, line=0.5)
  mtext(c("Homo Ref","Het","Homo Alt"), side=1,line=0,at=1:3)
  mtext(paste0("N=", bp$n), side=1,line=1,at=1:3)
  
  df$SNP=ifelse(as.numeric(as.character(df$SNP))==0,0,1)
  df$SNP=factor(df$SNP, levels=0:1)
  p0=signif(t.test(expression~SNP, df)$p.value,3)
  bp=boxplot(expression~SNP, data=df, ylab="Normalized Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
  stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE)
  title(main=paste0("dominant effect (pvalue=",p0,")"), cex.main=0.8, line=0.5)
  mtext(c("w/o allele","w/ allele"), side=1,line=0,at=1:2)
  mtext(paste0("N=", bp$n), side=1,line=1,at=1:2)
  
  mtext(paste("cis-eQTL for",g_symbol,"and",s2), outer = TRUE, cex = 1.5)        
}
dev.off() 