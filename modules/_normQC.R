###########################################
# Rscript to generate the RLE and NUSe plot to detect outlier
# Usage: Rscript $PATH/_normQC.R genes.fpkm.allSamples.uniq.xls genes.fpkm.allSamples.uniq.QC.pdf
# Reference: https://www.genevestigator.com/userdocs/manual/qc.html
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Aug-8
###########################################

args<-commandArgs(TRUE)

FPKMfile=args[1]  # either filename or stdin
outputfile=args[2]

fpkm=read.table(file(FPKMfile), header=T);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; fpkm=fpkm[,grep("FPKM",colnames(fpkm))]

# RLE: For each gene and each sample, ratios are calculated between the expression of a gene and the median expression of this gene across all samples of the experiment. For each sample, these relative expression values are displayed as a box plot. Since it is assumed that in most experiments only relatively few genes are differentially expressed, the boxes should be similar in range and be centered close to 0.

#Two effects may characterize arrays with lower quality: 1) the spread is greater than that of other arrays from this experiment, and 2) the box is not centered near 0.

pdf(outputfile, width=15, height=5)
fpkm=log10(fpkm+0.01)  # so row value of 0 will be -2 in the transformed value
fpkm=fpkm/apply(fpkm, 1, median)
colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))
require(reshape2)
fpkm=melt(cbind(ID=rownames(fpkm), fpkm), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(fpkm, reorder(Sample, FPKM, IQR))  # sort by IQR
boxplot(FPKM ~ bymedian, data=fpkm, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE (Relative Log Expression) plot", xlab="Sample ID", )

#The other graphical representation (NUSE) represents normalized standard error (SE) estimates from the PLM fit. The SE estimates are normalized such that for each probe set, the median standard error across all arrays is equal to 1. A box plot of NUSE values is drawn for each array. On the NUSE plot, arrays with lower quality will have boxes that are centered higher and/or have a larger spread than the other good quality arrays from the same experiment. Typically, boxes centered above 1.1 represent arrays that have quality problems which are often detected in several of the other QC measures presented in this chapter.

# TODO:

dev.off()