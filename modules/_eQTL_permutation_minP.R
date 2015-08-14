# Rscript to perform eQTL permutation with shuffled genotype labels
# Usage: Rscript _eQTL_permutation_minP.R 1 ~/neurogen/rnaseq_PD/results/eQTL/HCILBSNDA89samples/data.RData ~/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls.postPEER.xls
require(MatrixEQTL)

args<-commandArgs(TRUE)
i=args[1]
Rdata=args[2] # "~/neurogen/rnaseq_PD/results/eQTL/HCILBSNDA89samples/data.RData"
expr=args[3] # '~/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls.postPEER.xls'

load(Rdata)
residuals=read.table(expr, header=T)
genes = SlicedData$new();
genes$CreateFromMatrix(as.matrix(residuals))

message(paste(" -- performing the",i,"permutation  ..."))
snps_shuffled=snps$Clone();
snps_shuffled$ColumnSubsample(sample.int(ncol(snps_shuffled))) # scramble the sample lable of SNP

me2 = Matrix_eQTL_main(
  snps = snps_shuffled,
  gene = genes,
  cvrt = SlicedData$new(), 
  output_file_name = "",
  pvOutputThreshold = 0,  # no tran
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = FALSE,
  output_file_name.cis = paste0("final.cis.eQTL.permutation.K",i,".xls"),
  pvOutputThreshold.cis = 1e-5,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = 1e6,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = TRUE);

mpg=me2$cis$min.pv.gene
mpg=mpg[match(sort(names(mpg)), names(mpg))]
write.table(mpg,paste0("permutation",i,".txt"))