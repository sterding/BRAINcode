###############################################
## Rscript to run factor analysis using PEER
## Author: Xianjun Dong
## Date: 2014-Oct-27
## Version: 0.0
## Require: R v3.0.2, python etc.
## Usage: Rscript ~/neurogen/pipeline/RNAseq/modules/_factor_analysis.R
###############################################
require(reshape2)
require(peer)
require('MatrixEQTL') || install.package('MatrixEQTL', repo='http://cran.revolutionanalytics.com');

args<-commandArgs(TRUE)

expr_file=args[1]  # for example: expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls"
#gene_location_file_name = 
covariates_file_name=args[2]  # covariates_file_name="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"
snps_file_name=args[3]  # snps_file_name="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.txt"
snps_location_file_name=args[4]  # snps_location_file_name="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.SNP.ID"

# for test
expr_file="genes.fpkm.HCILB.uniq.80samples.xls";
covariates_file_name="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/covariances.12152014.80samples.tab";
snps_file_name="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.txt";
geneloc="genes.loci.txt";
snpsloc="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.SNP.ID"

message("# loading expression data...")
######################

# TODO: remove outlier and replicate
# TODO: change the sample ID to subject ID

expr = read.table(expr_file, header=T, check.names = F)
rownames(expr) = expr[,1]; expr = expr[, -1];
#expr=expr[,grep("FPKM", colnames(expr))]
#colnames(expr)=gsub("FPKM.","",colnames(expr))


message("# filtering out lowly expressed genes...")
######################
# remove genes with 0 in >=90% of samples
expr=expr[rowMeans(expr==0)<0.9, ]
# logorithm
expr=log10(expr+0.01)  # so row value of 0 will be -2 in the transformed value


message("loading SNP data...")
######################

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 5,000 rows
snps$LoadFile(snps_file_name);

snps$ColumnSubsample(colnames(snps) %in% colnames(expr))

message("# filter out SNPs with MAF<=0.05 ...");
# Note: here Minor allele is the minor one for the samples, not necessary the same one as the population)
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
cat('SNPs before filtering:',nrow(snps), "\n")  # 6109238
snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps), "\n")  # 6053947

## check
#df=read.table(snps_file_name, header=T, sep="\t", na.strings = "NA", check.names = F)
#rownames(df)=df[,1]; df=df[,-1]; 

message("# Load covariates ...")
######################

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

cvrt_snp=cvrt$Clone()
cvrt_snp$ColumnSubsample(grep("Genotype_batch",colnames(cvrt_snp)))
# covariates for genotyping can be used for eQTL. In our case, we discard it since we only have two replicates and the consistence is very high

cvrt$ColumnSubsample(grep("Genotype_batch",colnames(cvrt), invert =T))
cvrt$ColumnSubsample(grep("Diagnosis",colnames(cvrt), invert =T))

## TODO: change batch into categorical variables. e..g batch number has to be encoded as indicators, with a different binary variable for each batch, having a value of 1 if the individual was in the batch and a value of 0 otherwise. 

## two marriage strategies of PEER and Matrix-eQTL
######################
##1: include PEER factors as covariates
#- run PEER without including any covariates, 
#- get factors from PEER, e.g. factor = PEER_getX(model) 
#- use the factor as covariates for Matrix-eQTL
#
##2: regress them out beforehand
#- run PEER with covariates, e.g. PEER_setCovariates(model, as.matrix(covs))
#- get residuals from PEER, e.g. residuals = PEER_getResiduals(model), 
#- use the residual as expression input for Matrix-eQTL, but no more covariates included. It’s still OK to include covariates specific for genotypes, I guess.

    
# Method: http://geuvadiswiki.crg.es/index.php/Basic_methodology
#c. for each group, normalize with PEER, adding mean
#    c1. use subset (e.g. chr20, or chr20- 22) using K=0,1,3,5,7,10,13,15,20 for each dataset
#    c2. run eQTL and number of genes for each K.
#    c3. get the optimal K = K(with most number of eQTL genes)
#    c4. run PEER on 20,000 exons to get covairantes for the final normalization
#    c5. final PEER normalization using all dataset, residual + mean as final quantification
#d. transform the final quantification to standard normal distribution (by ?)
#e. eQTL using Matrix-eQTL: linear regression of quantification ~ genotypes + genotype_covariates


message("# loading SNP and gene position files...")
######################
snpspos = read.table(snpsloc, header = TRUE, stringsAsFactors = FALSE);
snpspos=snpspos[,1:3]
genepos = read.table(geneloc, header = TRUE, stringsAsFactors = FALSE);

message("# Run PEER on subset of genes to get best K ...")
######################
# random extract 5000 genes for test
expr_subset = expr[sample.int(nrow(expr), 5000),]  

nCiseqtl = c(); K=c(0,1,3,5,7,10,13,15,20);

for(k in K){
    model = PEER()
    PEER_setNk(model,k)
    PEER_setPhenoMean(model,as.matrix(t(expr_subset)))
    PEER_setAdd_mean(model, TRUE)
    PEER_setCovariates(model, t(as.matrix(cvrt)))
    PEER_getNk(model)
    PEER_update(model)
    residuals = t(PEER_getResiduals(model))
    rownames(residuals) = rownames(expr_subset)
    colnames(residuals) = colnames(expr_subset)
    
    # convert to SlideData format
    genes = SlicedData$new();
    genes$CreateFromMatrix(residuals);
    
    # convert to normal distribution
    for( sl in 1:length(genes) ) {
        mat = genes[[sl]];
        mat = t(apply(mat, 1, rank, ties.method = "average"));
        mat = qnorm(mat / (ncol(genes)+1));
        genes[[sl]] = mat;
    }
    rm(sl, mat);
    
    # run eQTL with residuals
    me = Matrix_eQTL_main(
        snps = snps,
        gene = genes,
        cvrt = SlicedData$new(),  # or cvrt_snp
        output_file_name = "",
        pvOutputThreshold = 0,  # no tran
        useModel = modelLINEAR, 
        errorCovariance = numeric(), 
        verbose = FALSE,
        output_file_name.cis = paste0("cis.eQTL.K",k,".txt"),
        pvOutputThreshold.cis = 1e-5,
        snpspos = snpspos, 
        genepos = genepos,
        cisDist = 1e6,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
    
    nCiseqtl = cbind(nCiseqtl, nrow(me$cis$eqtls))
    cat(k,":", nrow(me$cis$eqtls),"\n")
}

pdf("step1.bestK.pdf")
plot(K, nCiseqtl, type='b', xlab="PEER K", ylab="cis QTL genes", main="genes FPKM")
dev.off()

bestK = K[which.max(nCiseqtl)]

message("# now to get covariates ...")
######################
expr_subset = expr[sample.int(nrow(expr), 20000),]
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr_subset)))
PEER_setNk(model,bestK)
PEER_setAdd_mean(model, TRUE)
PEER_setCovariates(model, t(as.matrix(cvrt)))
PEER_update(model)
factors = PEER_getX(model)

message("# get residual...");
######################
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,bestK)
PEER_setCovariates(model, as.matrix(factors))
PEER_setAdd_mean(model, TRUE)
PEER_update(model)
residuals = t(PEER_getResiduals(model))
rownames(residuals) = rownames(expr)
colnames(residuals) = colnames(expr)

write.table(t(residuals), file = paste(expr_file, "peerResiduals.tab",sep="."), col.names=colnames(expr), row.names =F, sep = "\t", quote =F)

genes = SlicedData$new();
genes$CreateFromMatrix(residuals);

message("# convert to normal distribution")
######################

for( sl in 1:length(genes) ) {
    mat = genes[[sl]];
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(genes)+1));
    genes[[sl]] = mat;
}
rm(sl, mat);

message("# run eQTL with final quantifications")
######################

me = Matrix_eQTL_main(
    snps = snps,
    gene = genes,
    cvrt = SlicedData$new(), # or cvrt_snp,
    output_file_name = "final.trans.eQTL.txt",
    pvOutputThreshold = 1e-6,
    useModel = modelLINEAR, 
    errorCovariance = numeric(), , 
    verbose = TRUE,
    output_file_name.cis = "final.cis.eQTL.txt",
    pvOutputThreshold.cis = 1e-5,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

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