#!/usr/bin/env Rscript 
## Rscript to calculate the consistence of TNE expression across the samples
## Usage: Rscript /data/neurogen/pipeline/RNAseq/src/TNE_caller.consistency.R $SAMPLE_GROUP $list_bw_file 
## Author: Xianjun Dong
## Date: Nov-16-2017

args<-commandArgs(TRUE)

SAMPLE_GROUP = args[1]
list_bw_file = args[2]

list_bw = read.table(list_bw_file, header = F, col.names = c('sampleID','bigwigFile'), stringsAsFactors = F)

EXP=data.frame(); PV=data.frame(); 
for(bw_file in list_bw$bigwigFile){
  print(bw_file)
  # read background
  df=read.table(paste0(bw_file,".", SAMPLE_GROUP, ".rdbg"), header=F)[,2] # mean RPM (mean0 from bigWigAverageOverBed)
  Fn=ecdf(df)
  
  # read expression
  expression=read.table(paste0(bw_file,".", SAMPLE_GROUP, ".eRNA.meanRPM"), header=F)
  pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));

  # merge
  if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; }
  else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); }
}

colnames(EXP)=c("locus",list_bw$sampleID); colnames(PV)=c("locus",list_bw$sampleID); 
write.table(EXP, "eRNA.tmp5.meanRPM.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(PV,  "eRNA.tmp5.pvalues.xls", col.names=T, row.names=F, sep="\t", quote=F)

## binomial test for the significant HTNE (p<0.05)
N=nrow(list_bw)
binomial.pvalues = sapply(rowSums(PV[,-1]<=0.05), function(x) binom.test(x,N,0.05,'greater')$p.value)
# using the Holm-Bonferroni method (AKA step-down Bonferroni)  to correct for multiple test
p.adjusted = cbind(binomial.pvalues=binomial.pvalues, p.adjusted.HB = p.adjust(binomial.pvalues, method = "holm"), p.adjusted.bonferroni=p.adjust(binomial.pvalues, method = "bonferroni"),p.adjusted.FDR=p.adjust(binomial.pvalues, method = "fdr"))
rownames(p.adjusted) = PV[,1]

write.table(p.adjusted,  "eRNA.tmp5.pvalues.adjusted.xls", col.names=NA, row.names=T, sep="\t", quote=F)