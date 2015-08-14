###############################################
## Rscript to run factor analysis using PEER
## Author: Xianjun Dong
## Date: 2014-Dec-27
## Version: 0.0
## Require: R v3.0.2, python etc. e.g. module unload R/3.1.0; module load R/3.0.2
## Usage: bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_PEER_eQTL.R
###############################################
require(reshape2)
require(peer)
require(MatrixEQTL)

args<-commandArgs(TRUE)

expr_file=args[1]  # for example: expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls"
covs_file=args[2]  # covs_file="/PHShome/xd010/neurogen/rnaseq_PD/rawfiles/covariances.tab"
snps_file=args[3]  # snps_file="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.txt"
geneloc=args[4]
snpsloc=args[5]  # snps_location_file_name="/data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.SNP.ID"

covs_file="~/neurogen/rnaseq_PD/results/eQTL/HCILBSNDA89samples/covariance_tablel_109SNDAsubjects_v08072015.tab";
snps_file="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt";  # 93 unique subjects
snpsloc="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID"

# for gene
expr_file="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls";  # matrix of 57816 x 118 (HCILB for all cell types)
geneloc="~/neurogen/rnaseq_PD/results/merged/genes.loci.txt";

outliers=c('ILB_BN10-90_SNDA_4_rep1')

# for eRNAs
#expr_file="eRNA.RPKM.xls";  
#geneloc="eRNA.loci.txt";

if(file.exists("data.RData")) load("data.RData") else{
    
message("## loading expression data...")
######################
expr = read.table(expr_file, header=T, check.names = F)  # GxN where G is number of genes and N is number of samples
rownames(expr) = expr[,1]; expr = expr[, -1];
expr=expr[,grep("FPKM", colnames(expr))]
colnames(expr) = gsub("FPKM.","",colnames(expr))

message(" # remove outlier and replicate...")
######################
covs=read.table(covs_file, header=T)
covs=subset(covs, Diagnosis<2)  # only HC and ILB
covs=subset(covs, !(RNAseq_ID %in% outliers))  # remove outliers
expr=subset(expr, select = as.character(covs$RNAseq_ID))

# change the sample ID to subject ID
colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr))

message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))

message(" # filtering expression data...")
######################
# remove genes with 0 in >=90% of samples
# expr=expr[rowMeans(expr==0)<0.9, ]
# GTEx: Filter on >=10 individuals having >0.1 RPKM.
expr=expr[rowSums(expr>0.05)>=10,]  # 57816 --> 36556 remained
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

message("## loading SNP data...")
######################

snps = SlicedData$new();  # # GxN matrix
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 5,000 rows
snps$LoadFile(snps_file);

# extract the samples with both RNAseq and genotyped, and reorder both snps and expression
snps$ColumnSubsample(match(intersect(colnames(snps), colnames(expr)), colnames(snps)))
expr=expr[,intersect(colnames(snps), colnames(expr))]

message("# filter out SNPs with MAF<=0.05 ...");
# Note: here Minor allele is the minor one for the samples, not necessary the same one as the population)
maf.list = vector('list', length(snps))
na.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
  
  na.list[[sl]] = is.na(rowMeans(slice));
}
maf = unlist(maf.list)
na  = unlist(na.list)
cat('SNPs before filtering:',nrow(snps), "\n")  # 6109238
snps$RowReorder(!na & maf>0.05);  # remove rows including NA
cat('SNPs after filtering:',nrow(snps), "\n")  # 4320519

rm(maf, na, maf.list, na.list)

message("# Load covariates ...")
######################

cvrt = SlicedData$new();
rownames(covs)=covs$Subject_ID; covs=subset(covs, select = c(RNAseq_batch1, RNAseq_batch2, RNAseq_batch3,RNAseq_batch4, RNAseq_batch5, RNAseq_batch6, readsLength_75nt, RIN, Sex, Age, PMI))
cvrt$CreateFromMatrix(as.matrix(covs));

cvrt_snp=cvrt$Clone()
cvrt_snp$ColumnSubsample(grep("Genotype_batch",colnames(cvrt_snp)))
# covariates for genotyping can be used for eQTL. In our case, we discard it since we only have two replicates and the consistence is very high

#cvrt$ColumnSubsample(grep("Genotype_batch",colnames(cvrt), invert =T))
#cvrt$ColumnSubsample(grep("Diagnosis",colnames(cvrt), invert =T))

# sort the samples 
cvrt$RowReorder(match(intersect(colnames(snps),rownames(cvrt)), rownames(cvrt)))

message("# loading SNP and gene position files...")
######################
snpspos = read.table(snpsloc, header = TRUE, stringsAsFactors = FALSE);
snpspos=snpspos[,1:3]
genepos = read.table(geneloc, header = TRUE, stringsAsFactors = FALSE);

save.image("data.RData")
}

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
#- use the residual as expression input for Matrix-eQTL, but no more covariates included. It's still OK to include covariates specific for genotypes (I think).

# GEUVADIS used the #2 strategy, as below:
    
# Method: http://geuvadiswiki.crg.es/index.php/Basic_methodology#Normalization_of_the_phenotype_data
#c. for each group, normalize with PEER, adding mean
#    c1. use subset (e.g. chr20, or chr20- 22) using K=0,1,3,5,7,10,13,15,20 for each dataset
#    c2. run eQTL and number of genes for each K.
#    c3. get the optimal K = K(with most number of eQTL genes)
#    c4. run PEER on 20,000 exons to get covairantes for the final normalization
#    c5. final PEER normalization using all dataset, residual + mean as final quantification
#d. transform the final quantification to standard normal distribution (by ?)
#e. eQTL using Matrix-eQTL: linear regression of quantification ~ genotypes + genotype_covariates

message("# step1: getting the best K ...")
######################
# use genes on chr20-22
expr_subset = expr[rownames(expr) %in% subset(genepos, chr=="chr20"|chr=="chr21"|chr=="chr22")[,1], ]  # 2087 genes selected

n_eGene=c(); K=c(1,2,4,5,7,10,15,20);

for(k in K){
    model = PEER()
    PEER_setNk(model,k)
    PEER_setPhenoMean(model,as.matrix(t(expr_subset)))  # PEER ask NxG matrix, where N=samples and G=genes
    #PEER_setAdd_mean(model, TRUE)
    #PEER_setCovariates(model, as.matrix(cvrt))  # PEER ask NxC matrix, where N=samples and C=covariates
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)
    #PEER_getX(model)
    residuals = t(PEER_getResiduals(model))  # convert to GxN
    rownames(residuals) = rownames(expr_subset)
    colnames(residuals) = colnames(expr_subset)

#     pdf(file=paste("diagnostics_peer_", k, ".exprchr20_22.pdf", sep=""), paper="usr", width = 0, height = 0)
#     PEER_plotModel(model)
#     mtext(paste0("PEER factor = ", k), side=3, line=1)
#     dev.off()
    
    # add mean
    #residuals = residuals + apply(expr_subset, 1, mean)
    
    # convert to SlideData format
    genes = SlicedData$new();
    genes$CreateFromMatrix(residuals);
    
    rm(residuals, model)
    
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
        output_file_name.cis = paste0("step1.cis.eQTL.K",k,".xls"),
        pvOutputThreshold.cis = 1e-5,
        snpspos = snpspos, 
        genepos = genepos,
        cisDist = 1e6,
        pvalue.hist = 'qqplot',
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = TRUE);
    
    pdf(file=paste("diagnostics_matrixeqlt_", k, ".exprchr20_22.pdf", sep=""), paper="usr", width = 0, height = 0)
    plot(me, pch = 16, cex = 0.7);
    dev.off()
    
    n_eGene = cbind(n_eGene, sum(me$cis$min.pv.gene<1e-5))  
    cat(k,":", sum(me$cis$min.pv.gene<1e-5),"\n")
    
    rm(me)
}

pdf("step1.bestK.pdf")
plot(K, n_eGene, type='b', xlab="PEER K", ylab="# of cis-eGenes", main="pilot PEER + Matrix-eQTL")
dev.off()

bestK = K[which.max(n_eGene)]

message(paste0("\t bestK = ",bestK))
message(paste0("\t number of eGenes = ",max(n_eGene)))

message("# step2: getting covariates ...")
######################
#bestK=4
expr_subset = expr[sample.int(nrow(expr), 20000),]
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr_subset)))
PEER_setNk(model,bestK)
#PEER_setAdd_mean(model, TRUE)
#PEER_setCovariates(model, as.matrix(cvrt)) # include known cvrt as above
PEER_update(model)
factors = PEER_getX(model)  # now, the getX() factors should include mean + known covariates + learnt hidden factors.

cor(factors, covs)

rm(expr_subset, model)

message("# step3: getting residuals...");
######################
model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setNk(model,0)  #since we use above learnt factors as ?known? covariates, we don?t need to set number of hidden factors any more.
PEER_setAdd_mean(model, FALSE)  # mean is already included in the above getX(), so no need to include again.
PEER_setCovariates(model, as.matrix(factors))
PEER_update(model)
residuals = t(PEER_getResiduals(model))  # transfer to GxN matrix
rownames(residuals) = rownames(expr)
colnames(residuals) = colnames(expr)

pdf(file="diagnostics_peer.all.pdf", paper="usr", width = 0, height = 0)
PEER_plotModel(model)
mtext(paste0("PEER factor = ", k), side=3, line=1)
dev.off()

message("# add mean to residuals")
######################
#residuals = residuals + apply(expr, 1, mean)

#message("# convert to normal distribution")
#residuals = t(apply(residuals, 1, rank, ties.method = "average"));
#residuals = qnorm(residuals / (ncol(residuals)+1));

message("# run RLE on PEER normalized quantification data ...")
######################

## RLE before and after peer

pdf("RLE.plot.pdf", width=10, height=5)
res=data.frame(expr)
rle1=res/apply(res, 1, median)

res=data.frame(residuals)
rle2=res/apply(res, 1, median)

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

message("# save final quantification data into file")
######################
write.table(format(residuals, digits=4,nsmall=4), file = paste(expr_file, "postPEER.xls",sep="."), sep="\t", col.names = NA, quote=F,row.names = TRUE)

residuals=as.matrix(read.table(paste(expr_file, "postPEER.xls",sep="."), header=T))


message("# step4: run eQTL with final quantifications")
######################
genes = SlicedData$new();
genes$CreateFromMatrix(residuals);

me = Matrix_eQTL_main(
    snps = snps,
    gene = genes,
    cvrt = SlicedData$new(),
    output_file_name = "", #final.trans.eQTL.xls",  # we mute the trans-eQTL analysis at the moment
    pvOutputThreshold = 0, #1e-8,
    useModel = modelLINEAR, 
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = "final.cis.eQTL.xls",
    pvOutputThreshold.cis = 1e-6,  # a rough cutoff here in order to refine later
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = 1e6,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);

write.table(me$cis$min.pv.gene, file='final.cis.min.pv.gene.txt')

pdf(file="diagnostics_matrixeqlt.all.pdf", paper="usr", width = 0, height = 0)
plot(me, pch = 16, cex = 0.7);
dev.off()

message(paste("eQTL analysis is done and found", me$cis$neqtls,"cis SNP-gene pairs with significance p<1e-5."))

message("# step5: perform permutations for eGene FDR")
######################
snps_shuffled=snps$Clone();
min.pv.gene=c();
for(i in 1:1000)
{
  message(paste(" -- performing the",i,"permutation  ..."))
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
  
  if(length(min.pv.gene)==0) {
    min.pv.gene=mpg
  } else {
    min.pv.gene=cbind(min.pv.gene, mpg)
  }
}

# now min.pv.gene is a matrix for min(p) of all genes (per row) of 1000 permutations (per column)
observedP=me$cis$min.pv.gene

## alternatively, run _eQTL_permutation_minP.R parallelly and then merge together
setwd("~/neurogen/rnaseq_PD/results/eQTL/HCILBSNDA89samples/")
observedP=read.table('final.cis.min.pv.gene.txt')
n=sort(rownames(observedP))
observedP=observedP[match(sort(rownames(observedP)), rownames(observedP)),]
names(observedP)=n
min.pv.gene=read.table(paste0("permutation",1,".txt"))
for(i in 2:1000){
  message(i)
  min.pv.gene=cbind(min.pv.gene,read.table(paste0("permutation",i,".txt")))
}

#For each gene, adjusted P value is set to the rank of observed in the permutation list divided by X
adjustedP=apply(min.pv.gene < observedP, 1, mean)  # get the percenage of observedP>min(p) among the 1000 permutations

#Use adjusted P values as input to qvalue to estimate FDR for a given adjust P value cutoff.
require('qvalue') # source("http://bioconductor.org/biocLite.R"); biocLite("qvalue")
q=qvalue(p = adjustedP)$qvalue
write.table(q[q<=0.05],"final.cis.eGene.qvalue.cutoff0.05.txt")

## now, we can get a list of all significantly associated SNP-gene pairs
#get the p-value correspinding to the 0.05 tail of the distribution of the minimal p-values from 1000 permutations 
#since the distribution is asymmetric, I will simply use the 5 percentile of sorted values 
cutoff = apply(min.pv.gene, 1, quantile, probs=0.05)

genesnp=read.table("final.cis.eQTL.xls", header=T, stringsAsFactors=F)
genesnp = genesnp[with(genesnp, order(FDR)), ]
genesnp = subset(genesnp, FDR<0.05)

## if any, add gene symbol etc. annotation info to the eQTL output
annotation = read.table("/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header=F)
colnames(annotation) = c("chr","start", "end", "gene_symbol", "gene_ensID", "gene_type","strand")
genesnp = cbind(genesnp,annotation[match(genesnp$gene, annotation$gene_ensID),])

eGenes=aggregate(SNP~gene_symbol+gene_ensID+gene_type, data=data.frame(genesnp), FUN=length)
write.table(eGenes[order(-eGenes$SNP),], file="eGenes.FDR.05.sortbyeQTL.xls", sep="\t", col.names = T, quote=F,row.names = F)

write.table(genesnp, file="final.cis.eQTL.FDR.05.xls", sep="\t", col.names = T, quote=F,row.names = F)

message("# making eQTL plot ...")
######################
# load data
load("data.RData")
genesnp = read.table("final.cis.eQTL.FDR.05.xls", header=T, stringsAsFactors =F)
expr_file="/PHShome/xd010/neurogen/rnaseq_PD/results/eQTL/genes80samples/genes.fpkm.HCILB.uniq.80samples.xls";  # matrix of 57816 x 80
residuals = read.table(paste(expr_file, "postPEER.xls",sep="."))
require(MatrixEQTL)
genes = SlicedData$new();
genes$CreateFromMatrix(as.matrix(residuals))

# one file per gene
for(g in unique(genesnp$gene))
{
    pdf(paste(expr_file, g, "eQTLplot.pdf", sep="."), width=8, height=8)
    par(mfrow=c(1,2), mar=c(4,4,2,2), oma = c(0, 0, 2, 0))
    genesnp0=subset(genesnp, gene==g)
    genesnp0 = genesnp0[with(genesnp0, order(FDR)), ]
    message(g);
    for(i in 1:nrow(genesnp0))
    {
        s=as.character(genesnp0$SNP[i]);
        s2=sub("([a-zA-Z].*):(.*):(.*):(.*)","\\1", s) 
        g=genesnp0$gene[i];
        g_symbol=ifelse(is.null(genesnp0$gene_symbol[i]), genesnp0$gene[i], genesnp0$gene_symbol[i]);
        p=signif(genesnp0$p.value[i], 3);
        print(paste(i, s, s2, g_symbol, p))
        
        df=data.frame(expression=as.numeric(genes$FindRow(g)$row), SNP=as.numeric(snps$FindRow(s)$row))
        df$SNP=factor(df$SNP, levels=0:2)
        bp=boxplot(expression~SNP, data=df, ylab="Rank Normalized Gene Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
        stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE) 
        title(main=paste0("additive effect (pvalue=", p,")"), cex.main=0.8, line=0.5)
        mtext(c("Homo Ref","Het","Homo Alt"), side=1,line=0,at=1:3)
        mtext(paste0("N=", bp$n), side=1,line=1,at=1:3)
    
        df$SNP=ifelse(as.numeric(as.character(df$SNP))==0,0,1)
        df$SNP=factor(df$SNP, levels=0:1)
        p0=signif(t.test(expression~SNP, df)$p.value,3)
        bp=boxplot(expression~SNP, data=df, ylab="Rank Normalized Gene Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
        stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE)
        title(main=paste0("dominant effect (pvalue=",p0,")"), cex.main=0.8, line=0.5)
        mtext(c("w/o allele","w/ allele"), side=1,line=0,at=1:2)
        mtext(paste0("N=", bp$n), side=1,line=1,at=1:2)
        
        mtext(paste("cis-eQTL for",g_symbol,"and",s2), outer = TRUE, cex = 1.5)        
    }
    dev.off()    
}

# one file for all genes
pdf(paste(expr_file, "eQTLplot.pdf", sep="."), width=8, height=8)
par(mfrow=c(1,2), mar=c(4,4,2,2), oma = c(0, 0, 2, 0))
for(i in 1:1000){
    s=as.character(genesnp$SNP[i]);
    #clean up SNP name (for those with rs6660464:72826949:T:C)
    s2=sub("([a-zA-Z].*):(.*):(.*):(.*)","\\1", s) 
    g=genesnp$gene[i];
    g_symbol=ifelse(is.null(genesnp$gene_symbol[i]), genesnp$gene[i], genesnp$gene_symbol[i]);
    p=signif(genesnp$p.value[i], 3);
    print(paste(i, s, s2, g_symbol, p))
    
    df=data.frame(expression=as.numeric(genes$FindRow(g)$row), SNP=as.numeric(snps$FindRow(s)$row))
    df$SNP=factor(df$SNP, levels=0:2)
    bp=boxplot(expression~SNP, data=df, ylab="Rank Normalized Gene Expression log10(RPKM)", xaxt='n', main="",  col='lightgreen', outpch=NA)
    stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE) 
    title(main=paste0("additive effect (pvalue=", p,")"), cex.main=0.8, line=0.5)
    mtext(c("Homo Ref","Het","Homo Alt"), side=1,line=0,at=1:3)
    mtext(paste0("N=", bp$n), side=1,line=1,at=1:3)

    df$SNP=ifelse(as.numeric(as.character(df$SNP))==0,0,1)
    df$SNP=factor(df$SNP, levels=0:1)
    p0=signif(t.test(expression~SNP, df)$p.value,3)
    bp=boxplot(expression~SNP, data=df, ylab="Rank Normalized Gene Expression log10(RPKM)", xaxt='n', main="", col='lightblue', outpch=NA)
    stripchart(expression~SNP, data=df, vertical = TRUE, method = "jitter", pch = 1, col = "darkred", cex=0.6, add = TRUE)
    title(main=paste0("dominant effect (pvalue=",p0,")"), cex.main=0.8, line=0.5)
    mtext(c("w/o allele","w/ allele"), side=1,line=0,at=1:2)
    mtext(paste0("N=", bp$n), side=1,line=1,at=1:2)
    
    mtext(paste("cis-eQTL for",g_symbol,"and",s2), outer = TRUE, cex = 1.5)
}
dev.off()