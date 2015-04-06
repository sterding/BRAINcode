###############################################
## Rscript to run PEER normalization
## Author: Xianjun Dong
## Date: 2015-02-26
## Version: 0.0
## Require: R v3.0.2, python etc.
## Usage: bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_PEER_normalize.R
###############################################


args<-commandArgs(TRUE)

expr_file=args[1]  # expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls"
covs_file=args[2]  # covs_file="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"

expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls"
covs_file="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"
bestK = 10 # according to eQTL test using SNDA data


if(file.exists("data.RData")) load("data.RData") else{
    
message("# loading expression data...")
######################

# TODO: remove outlier and replicate

expr = read.table(expr_file, header=T, check.names = F, stringsAsFactors = F)  # GxN where G is number of genes and N is number of samples
rownames(expr) = expr[,1]; expr = expr[, -1];
expr=expr[,grep("FPKM", colnames(expr))]
# remove FPKM in the colname
colnames(expr)=gsub("FPKM.(.*)", "\\1", colnames(expr))

# change the sample ID to subject ID
#colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr))

message("# Load covariates ...")
######################

cvrt = read.table(covs_file, header=T, check.names = F, stringsAsFactors = F) # NxC matrix
cvrt = subset(cvrt, outlier==0, select=-c(outlier, subjectID, replicate))  # remove outlier

# Encode categorical variables (e.g., batch, sex, readsLength, condition, cellType, replicate) as indicators. E.g. with a different binary variable for each batch, having a value of 1 if the individual was in the batch and a value of 0 otherwise.
library(reshape2)
categorical_varaibles = c("cellType", "batch", "sex", "readsLength", "condition");
for(x in categorical_varaibles) {cvrt = cbind(cvrt, value=1); cvrt[,x]=paste0(x,cvrt[,x]); cvrt = dcast(cvrt, as.formula(paste0("... ~ ", x)), fill=0);}

# intersected samples in covariates table & RNAseq
rownames(cvrt) = cvrt[, 1]; cvrt = cvrt[, -1]; 
expr=expr[, intersect(rownames(cvrt), colnames(expr))]
cvrt=cvrt[intersect(rownames(cvrt), colnames(expr)), ]

# temp: remove condition
cvrt = cvrt[,grep("condition", colnames(cvrt), invert=T)]

message("# filtering expression data...")
######################

# logorithm
expr=log10(expr+1e-5)  # so row value of 0 will be -5 in the transformed value

# like GTEx: Filter on >=10 individuals having >0.05 RPKM.
# expr=expr[rowSums(expr>0.05)>=10,]  # 57816 --> 36556 remained
pdf("individual.hist.beforePEER.beforeFilter.pdf",paper="usr", width = 0, height = 0)
par(mfrow=c(5,5), mar=c(2,2,2,1))
sapply(colnames(expr),function(x) hist(expr[,x],main=x,breaks=100, cex.main=0.8))
replicate(ceiling(ncol(expr)/5)*5-ncol(expr),plot.new())
sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"),function(x) hist(apply(expr[,grep(x,colnames(expr))],1,median),main=x,breaks=100, cex.main=0.9))
dev.off()

## Filter criteria: at least 2 samples or 1/2 of samples, whichever is bigger, having expression >0 RPKM.
#index= apply(expr[,grep("MCPY",colnames(expr))]>-15,1,sum)>=2 | apply(expr[,grep("TCPY",colnames(expr))]>-15,1,sum)>=5 | apply(expr[,grep("SNDA",colnames(expr))]>-15,1,sum)>=44 
## at least 10 samples have expression > 10^-15.
#index= rowSums(expr > -15)>=10  # 57816 --> 41730

# TODO: change to require median > -4 in at least one cell type (see distribution for how to define the cutoff of -4). # 57816 --> 34194
index= apply(expr[,grep("SNDA",colnames(expr))],1,median) >= -4 | apply(expr[,grep("TCPY",colnames(expr))],1,median) >= -4 | apply(expr[,grep("MCPY",colnames(expr))],1,median) >= -4

# how to define cutoff -4
pdf("celltype.hist.median.cutoff.pdf", paper='usr', width=0, height=4)
par(mfrow=c(1,3))
lapply(c("SNDA","TCPY","MCPY"),function(x) hist(apply(expr[,grep(x, colnames(expr))], 1, median), border=NA, col='darkblue', xlab="median of log10(FPKM)", breaks=200, main=x))
lapply(c("SNDA","TCPY","MCPY"),function(x) {hist(apply(expr[,grep(x, colnames(expr))], 1, median), border=NA, col='darkblue', xlab="median of log10(FPKM)",breaks=200, xlim=c(-8,5), ylim=c(0,1600), main=x); abline(v=-4,col='red',lwd=3)})
lapply(c("SNDA","TCPY","MCPY"),function(x) hist(apply(expr[index,grep(x, colnames(expr))], 1, median), border=NA, col='darkblue', xlab="median of log10(FPKM)",breaks=200, main=x))
dev.off()

expr = expr[index, ]



pdf("individual.hist.beforePEER.afterFilter.pdf",paper="usr", width = 0, height = 0)
par(mfrow=c(5,5), mar=c(2,2,2,1))
sapply(colnames(expr),function(x) hist(expr[,x],main=x,breaks=100, cex.main=0.8))
replicate(ceiling(ncol(expr)/5)*5-ncol(expr),plot.new())
sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"),function(x) hist(apply(expr[,grep(x,colnames(expr))],1,median),main=x,breaks=100, cex.main=0.9))
dev.off()

# outlier correction: quantile normalization with order preserved. Now RPKM is changed to rank normalized gene expression.
m=apply(expr, 1, mean); sd=apply(expr, 1, sd)
expr1 = t(apply(expr, 1, rank, ties.method = "average"));
#expr = qnorm(expr / (ncol(expr)+1));  # to standard normalization
expr1 = qnorm(expr1 / (ncol(expr1)+1), mean=m, sd=sd)  # or, to preserve the mean and sd of each gene

pdf("individual.hist.beforePEER.afterFilter.afterQuanNorm.pdf",paper="usr", width = 0, height = 0)
par(mfrow=c(5,5), mar=c(2,2,2,1))
sapply(colnames(expr1),function(x) hist(expr1[,x],main=x,breaks=100, cex.main=0.8))
replicate(ceiling(ncol(expr1)/5)*5-ncol(expr1),plot.new())
sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"),function(x) hist(apply(expr1[,grep(x,colnames(expr1))],1,median),main=x,breaks=100, cex.main=0.9))
dev.off()

rm(m,sd)

save.image("data.RData")
}


message("# Run PEER to normalize the expression ...")
######################
require(peer)
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr1)))
PEER_setNk(model, 0)  # no hidden factors
#PEER_setAdd_mean(model, TRUE)  #
PEER_setCovariates(model, as.matrix(cvrt)) # include known cvrt as above
PEER_setNmax_iterations(model, 3000)
PEER_update(model)
X = PEER_getX(model)  # now, the getX() factors should include known covariates + mean
W = PEER_getW(model)

residuals = t(PEER_getResiduals(model))  # convert to GxN
rownames(residuals) = rownames(expr1)
colnames(residuals) = colnames(expr1)

# explained variance by cell type
Yp = W[, -c(4:6)] %*% t(X[, -c(4:6)]) ## since the 4:6 columns are for celltypes. 
res1 = expr1 -Yp  # residuals of cell type

## solution (A): add the overall mean to the residuals, which could potentially flat the difference between cell types, since the residulas are trival comparing to overall mean. 
res1a = res1 + apply(expr1, 1, mean)

# solution (B): add the celltype-specific mean to the residuals
x=apply(expr1[, grep("SNDA",colnames(expr1))], 1, mean)
res1b=res1[,grep("SNDA", colnames(expr1))] + x
x=apply(expr1[, grep("TCPY",colnames(expr1))], 1, mean)
res1b=cbind(res1b, res1[,grep("TCPY", colnames(expr1))] + x)
x=apply(expr1[, grep("MCPY",colnames(expr1))], 1, mean)
res1b=cbind(res1b, res1[,grep("MCPY", colnames(expr1))] + x)
colnames(res1b) = colnames(expr1)[c(grep("SNDA",colnames(expr1)), grep("TCPY",colnames(expr1)), grep("MCPY", colnames(expr1)))]

# explained variance by batch effect
Yp2 = W[, -c(7:11)] %*% t(X[, -c(7:11)])  
res2 = expr1 -Yp2  # residuals of batch effect

pdf("PEER.effects.on.batch.and.celltype.pdf", paper='usr', width=0)
par(mfrow=c(1,6), mar=c(5,2,3,1), tck=0.03, las=2)
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(expr[,grep(x,colnames(expr))], 1, mean)}),   ylim=c(-5,0), main="raw RPKM")
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(residuals[,grep(x,colnames(residuals))], 1, mean)}), ylim=c(-1,1), main="PEER residuals\n(adjust all)")
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(res1[,grep(x,colnames(res1))], 1, mean)}),  ylim=c(-1,1), main="PEER residuals\n(adjust all but cellType)")
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(res1a[,grep(x,colnames(res1a))], 1, mean)}),    ylim=c(-5,0), main="res + overallMean\n(adjust all but cellType)")
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(res1b[,grep(x,colnames(res1b))], 1, mean)}),  ylim=c(-5,0), main="res + cellMean\n(adjust all but cellType)")
boxplot(sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"), function(x) {apply(res2[,grep(x,colnames(res2))], 1, mean)}),  ylim=c(-1,1), main="PEER residuals\n(adjust all but batch)")
dev.off()

pdf("individual.hist.afterPEER.pdf",paper="usr", width = 0, height = 0)
par(mfrow=c(5,5), mar=c(2,2,2,1))
sapply(colnames(res1b),function(x) hist(res1b[,x],main=x,breaks=100, cex.main=0.8))
replicate(ceiling(ncol(res1b)/5)*5-ncol(res1b),plot.new())
sapply(c(paste0("SNDA_",1:5), "MCPY_4", "TCPY_5"),function(x) hist(apply(res1b[,grep(x,colnames(res1b))],1,median),main=x,breaks=100, cex.main=0.9))
dev.off()

message("# save final quantification data into file")
######################
write.table(format(res1b, digits=4,nsmall=4), file = paste(expr_file, "postPEER.addGroupmean.xls",sep="."), sep="\t", col.names = NA, quote=F,row.names = TRUE)
write.table(format(res1, digits=4,nsmall=4), file = paste(expr_file, "postPEER.residuals.xls",sep="."), sep="\t", col.names = NA, quote=F,row.names = TRUE)

# check batch effect 
expr0=expr[,grep("SNDA", colnames(expr))]
require(peer)
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr0)))
PEER_setNk(model, 0)  # no hidden factors
PEER_setAdd_mean(model, TRUE)  #
PEER_setCovariates(model, as.matrix(subset(cvrt, cellTypeSNDA==1, select=c(RIN, age, PMI, batch1, batch2, batch3, batch4, batch5))))
PEER_setNmax_iterations(model, 3000)
PEER_update(model)
X = PEER_getX(model)  # now, the getX() factors should include known covariates + mean
W = PEER_getW(model)
residuals = t(PEER_getResiduals(model))  # convert to GxN
rownames(residuals) = rownames(expr0)
colnames(residuals) = colnames(expr0)

res0 = residuals + apply(expr0, 1, mean)

res1 = expr0 - W %*% t(X)
rownames(res1) = rownames(expr0)
colnames(res1) = colnames(expr0)

res2 = expr0 - W[, -c(4:8)] %*% t(X[, -c(4:8)])
rownames(res2) = rownames(expr0)
colnames(res2) = colnames(expr0)

res3 = expr0 - W[, c(4:8)] %*% t(X[, c(4:8)])
rownames(res3) = rownames(expr0)
colnames(res3) = colnames(expr0)

par(mfrow=c(1,3))
boxplot(sapply(1:5, function(x) {apply(residuals[,grep(paste0("SNDA_",x),colnames(residuals))], 1, mean)}), ylim=c(-3,4), main="PEER_getResiduals")
boxplot(sapply(1:5, function(x) {apply(expr0[,grep(paste0("SNDA_",x),colnames(expr0))], 1, mean)}), ylim=c(-3,4), main="raw RPKM")
boxplot(sapply(1:5, function(x) {apply(res2[,grep(paste0("SNDA_",x),colnames(res2))], 1, mean)}), ylim=c(-3,4), main="expr0 - W[,-c(4:8)]*X[,-c(4:8)]")



## without add mean effect

expr0=expr[,grep("SNDA", colnames(expr))]
require(peer)
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr0)))
PEER_setNk(model, 0)  # no hidden factors
PEER_setAdd_mean(model, FALSE)  #
PEER_setCovariates(model, as.matrix(subset(cvrt, cellTypeSNDA==1, select=c(RIN, age, PMI, batch1, batch2, batch3, batch4, batch5))))
PEER_setNmax_iterations(model, 3000)
PEER_update(model)
X = PEER_getX(model)  # now, the getX() factors should include known covariates + mean
W = PEER_getW(model)
residuals = t(PEER_getResiduals(model))  # convert to GxN
rownames(residuals) = rownames(expr0)
colnames(residuals) = colnames(expr0)

res0 = residuals + apply(expr0, 1, mean)

res1 = expr0 - W %*% t(X)
rownames(res1) = rownames(expr0)
colnames(res1) = colnames(expr0)

res2 = expr0 - W[, -c(4:8)] %*% t(X[, -c(4:8)])
rownames(res2) = rownames(expr0)
colnames(res2) = colnames(expr0)

res3 = W[, c(4:8)] %*% t(X[, c(4:8)])
rownames(res3) = rownames(expr0)
colnames(res3) = colnames(expr0)

par(mfrow=c(1,5))
boxplot(sapply(1:5, function(x) {apply(residuals[,grep(paste0("SNDA_",x),colnames(residuals))], 1, mean)}), ylim=c(-3,4), main="PEER_getResiduals")
boxplot(sapply(1:5, function(x) {apply(res1[,grep(paste0("SNDA_",x),colnames(res1))], 1, mean)}), ylim=c(-3,4), main="expr0 - W*X")
boxplot(sapply(1:5, function(x) {apply(expr0[,grep(paste0("SNDA_",x),colnames(expr0))], 1, mean)}), ylim=c(-3,4), main="raw RPKM")
boxplot(sapply(1:5, function(x) {apply(res2[,grep(paste0("SNDA_",x),colnames(res2))], 1, mean)}), ylim=c(-3,4), main="expr0 - W[,-c(4:8)]*X[,-c(4:8)]")
boxplot(sapply(1:5, function(x) {apply(res3[,grep(paste0("SNDA_",x),colnames(res3))], 1, mean)}), ylim=c(-3,4), main="W[,4:8]*X[,4:8]")


message("# run RLE on PEER normalized quantification data ...")
######################

## RLE before and after peer

pdf("RLE.plot.pdf", width=10, height=5)
res=data.frame(expr)
rle1=res/apply(res, 1, median)

res=data.frame(residuals)
rle2=res/apply(res, 1, median)

require(reshape2)
rle=melt(cbind(ID=rownames(rle1), rle1), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-2.5,4.5), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot before PEER", xlab="", ylab="Relative log expression (RLE)")
abline(h=1, col='red',lwd=1)

## RLE after peer
rle=melt(cbind(ID=rownames(rle2), rle2), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-2.5,4.5), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot after PEER", xlab="", ylab="Relative log expression (RLE)")
abline(h=1, col='red',lwd=1)

dev.off()

# expression distribution

pdf("expression.hist.plot.pdf", width=8, height=8)
par(mfrow=c(2,1))
hist(apply(expr,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution before PEER")
hist(apply(residuals,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution after PEER")
dev.off()

