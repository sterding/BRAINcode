# Rscript to calculate fold-change and p-value for HTNEs cross sample groups
args<-commandArgs(TRUE)
exp=args[1]
ref=args[2]
tars=args[3:length(args)]

# exp="~/eRNAseq/HCILB_SNDA/eRNA.meanRPM.allSamples.xls"
# ref="~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA"
# tars=c("~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY", "~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron")

EXP=read.table(exp, header=T, check.names = F, stringsAsFactors =F)
rownames(EXP) = EXP[,1]; EXP = EXP[,-1];
head(EXP)

REF=read.table(ref, header=F, check.names = F, stringsAsFactors =F)[,1]
ref=EXP[,REF]+0.001

RESULT=c();

for(i in tars){
  TAR=read.table(i, header=F, check.names = F, stringsAsFactors =F)[,1]
  tar=EXP[,TAR]+0.001
  results = apply(cbind(ref, tar), 1, function(x) {
    log2fc=log2(mean(x[1:length(REF)]) / mean(x[(length(REF)+1):(length(REF)+length(TAR))]));
    pvalue=t.test(x[1:length(REF)], x[(length(REF)+1):(length(REF)+length(TAR))], alternative="two.sided", var.equal=FALSE)$p.value;
    return(c(log2fc, pvalue));
    })
  results = t(results)
  colnames(results) = paste0(sub(".*\\.","",i),".", c('log2fc','pvalue'))
  RESULT=cbind(RESULT, results)
}
rownames(RESULT) = rownames(ref)

write.table(RESULT, sub("xls","ttest.txt",exp), col.names=T, row.names=T, sep="\t", quote=F)
