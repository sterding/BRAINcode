###############################################
## Rscript to run factor analysis using PEER
## Author: Xianjun Dong
## Date: 2014-Oct-27
## Version: 0.0
## Require: R v3.0.2, python etc. 
###############################################
require(reshape2)
library(peer)

args<-commandArgs(TRUE)

expr_file=args[1]  # for example: expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls"
covs_file=args[2]  # covs_file="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"

# for test
# expr_file="genes.fpkm.HCILB.uniq.xls"; covs_file="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"

message("loading data...")
######################

expr = read.table(expr_file, header=T)
rownames(expr) = expr[,1]; expr = expr[, -1];
expr=expr[,grep("FPKM", colnames(expr))]
colnames(expr)=gsub("FPKM.","",colnames(expr))
# add rep1 for earlier batch
colnames(expr) = ifelse(grepl("rep", colnames(expr)), colnames(expr), paste(colnames(expr),"_rep1",sep=""))

message("filtering out lowly expressed genes...")
######################
# remove genes with 0 in >=90% of samples
expr=expr[rowMeans(expr==0)<0.9, ]
# logorithm
expr=log10(expr+0.01)  # so row value of 0 will be -2 in the transformed value

#c. for each group, normalize with PEER, adding mean
#    c1. use subset (??e.g. chr20, or chr20- 22??) using K=????0,1,,3,,57,10,13,15,20 for each dataset
#    c2. run eQTL and number of genes for each K.
#    c3. get the optimal K = K(with most number of eQTL genes)
#    c4. run PEER on 20,000 exons to get covairantes for the final normalization
#    c5. final PEER normalization using all dataset, residual + mean as final quantification
#d. transform the final quantification to standard normal distribution (by ?)
#e. eQTL using Matrix-eQTL: linear regression of quantification ~ genotypes + genotype_covariates


# PEER without covariate
######################

model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,10)  # say we want to infer K=10 hidden confounders,
PEER_getNk(model)
PEER_update(model)

# output
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
plot(precision)
PEER_plotModel(model)

write.table(t(residuals), file = paste(expr_file, "peerResiduals.tab",sep="."), col.names=colnames(expr), row.names =F, sep = "\t", quote =F)

## RLE after peer
######################
rle=as.data.frame(t(residuals))
colnames(rle)=colnames(expr); rownames(rle)=rownames(expr)
rle=rle/apply(rle, 1, median)
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot after PEER", xlab="")

# PEER with covariates
######################

# read covariates table
covs = read.table(covs_file, header=T)
rownames(covs) = covs[,1]; covs = covs[, -1];
# TOFINISH: melt to long table  
covs=melt(cbind(ID=rownames(covs), covs), variable.name = "Sample",value.name ="FPKM", id="ID")

PEER_setCovariates(model, as.matrix(covs))
PEER_update(model)

plotcor(PEER_getX(model)[,1], PEER_getCovariates(model)[,1])
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
plot(precision)
PEER_plotModel(model)

## RLE after peer
######################
rle=as.data.frame(t(residuals))
colnames(rle)=colnames(expr); rownames(rle)=rownames(expr)
rle=rle/apply(rle, 1, median)
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot after PEER", xlab="")

dev.off()