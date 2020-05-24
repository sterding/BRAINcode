###############################################
## Rscript to run factor analysis using SVA
## Author: Xianjun Dong
## Date: 2016-Mar-20
## Version: 0.0
## Usage: bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_SVA.normalized.R samplelist.txt
###############################################
require(reshape2)
require(RCurl)
require(caret)
# source("http://bioconductor.org/biocLite.R"); biocLite("sva")
library(sva)

args<-commandArgs(TRUE)
SAMPLE_GROUP=args[1]

samplelist=read.table(paste0("~/neurogen/rnaseq_PD/results/merged/samplelist.",SAMPLE_GROUP), header = F, stringsAsFactors = FALSE)[,1]
covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples

expr_file="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls";  

message("## loading expression data...")
######################
expr = read.table(expr_file, header=T, check.names = F)  # GxN where G is number of genes and N is number of samples
rownames(expr) = expr[,1]; expr = expr[, -1];
if(grepl("cufflinks", expr_file)) expr=expr[,grep("FPKM", colnames(expr))]
colnames(expr) = gsub("FPKM.","",colnames(expr))

message(" # remove outlier and replicate...")
######################
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
#covs=subset(covs, BRAINCODE.final.selection==1)
common=intersect(samplelist, covs$sampleName)
covs=covs[match(common, covs$sampleName), ]
expr=subset(expr, select = common)

# change the sample ID to subject ID
colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr))

# convert to factor
covs$batch = as.factor(covs$batch)
covs$sex = as.factor(covs$sex)
covs$readsLength = as.factor(covs$readsLength)

message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))

message(" # filtering expression data...")
######################
# remove genes with 0 in >=90% of samples
# expr=expr[rowMeans(expr==0)<0.9, ]
# GTEx: Filter on >=10 individuals having >0.1 RPKM.
expr=expr[rowSums(expr>0.1)>=10,]  # 51506 --> 32642 remained
message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))

message(" # transforming RPKM to rank normalized gene expression ...")
######################
# logorithm
expr=log10(expr+0.01)  # so row value of 0 will be -2 in the transformed value
# outlier correction: quantile normalization with order preserved. Now RPKM is changed to rank normalized gene expression.
m=apply(expr, 1, mean); sd=apply(expr, 1, sd)
expr = t(apply(expr, 1, rank, ties.method = "average"));
#expr = qnorm(expr / (ncol(expr)+1));  # to standard normalization
expr = qnorm(expr / (ncol(expr)+1), mean=m, sd=sd)  # or, to preserve the mean and sd of each gene

rm(m,sd)


message("# adjusting expression with covariates...")
######################

# sva adjusted
Mod = model.matrix(~batch+sex+RIN+age+PMI, data=covs)
#Mod0 = model.matrix(~1,data=covs)

## estimate the latent factors
#n.sv = num.sv(expr, Mod, method="leek")
#svaobj = sva(as.matrix(expr),Mod, Mod0, n.sv=n.sv)
svaobj = sva(as.matrix(expr),Mod)
fsvaobj = fsva(dbdat=as.matrix(expr),mod=Mod,sv=svaobj, newdat=as.matrix(expr))

residuals = fsvaobj$db

message("# run RLE on SVA normalized quantification data ...")
######################

## RLE before and after SVA

pdf(paste(SAMPLE_GROUP,"RLE.plot.pdf",sep="."), width=10, height=5)
res=data.frame(expr, check.names = F)
names(res) = covs[match(covs$subjectID, names(res)),'sampleName']
rle1=res-apply(res, 1, median)

res=data.frame(residuals, check.names = F)
names(res) = covs[match(covs$subjectID, names(res)),'sampleName']
rle2=res-apply(res, 1, median)

rle=melt(cbind(ID=rownames(rle1), rle1), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-max(abs(range(rle$FPKM))),max(abs(range(rle$FPKM)))), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot before adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)

## RLE after SVA
rle=melt(cbind(ID=rownames(rle2), rle2), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-max(abs(range(rle$FPKM))),max(abs(range(rle$FPKM)))), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot after adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)

dev.off()

# expression distribution

pdf(paste(SAMPLE_GROUP,"expression.hist.plot.pdf",sep="."), width=8, height=8)
par(mfrow=c(2,1))
hist(apply(expr,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution before adjustment")
hist(apply(residuals,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution after adjustment")
dev.off()

message("# save final quantification data into file")
######################
write.table(format(residuals, digits=4,nsmall=4), file = paste(SAMPLE_GROUP,"expression.postSVA.xls",sep="."), sep="\t", col.names = NA, quote=F,row.names = TRUE)
