#!/usr/bin/Rscript
message("################################################################################################")
message("## Rscript to run PEER factor analyses with GGtools to enrich eQTL hits                        ##")
message("##  Author: Dr. Alex Tzuu-Wang Chang                                                            ##")
message("##    Date: 2015-07-30                                                                           ##")
message("## Version: 1.2 (added the codes for counting eGenes and more run-time information output)        ##")
message("## Require: R v3.2.1, python etc.                                                                  ##")
message("##                                                                                                  ##") 
message("## When submit jobs at the outside of Channing:                                                      ##")
message("## (1) Login to capecod first (2) qsub -l lx yourScript.sh                                            ##")
message("## Switch to LINUX node: ssh -l retwc -o StrictHostKeyChecking=no nantucket.bwh.harvard.edu            ##")
message("##                                                                                                      ##")
message("## Run R in Channing:[ qs | (qrsh -l lx6,large) | (qrsh -l lx) ] then R; qrsh -l rstudio; R --vanilla;  ##")
message("## /udd/stvjc/VM/R-devel-dist/bin/R; ~stvjc/bin/Rdevel --vanilla;                                       ##")
message("## Run R script in Channing:                                                                           ##")
message("## nohup R --no-save < myRscript.R > myRscriptR.out &                                                 ##")
message("## R --vanilla < /udd/retwc/R/test.R > /udd/retwc/R/testR.out ; qsub -l lx ./runMyRscript.sh        ##")
message("##                                                                                                  ##")
message("## After qs or qrsh the R executable file is located at the following path:                        ##")
message("## /local/bin/R -> /app/R-3.2.0@i86-rhel6.0/bin/R                                                  ##")
message("## Rscript -> /app/R-3.2.0@i86-rhel6.0/bin/Rscript                                                  ##")
message("## Check Job status: qds retwc; qstat -r; qstat -f; qstat -f | grep retwc | less; qstat -ls; qhost;  ##")
message("##                                                                                                    ##")
message("## Usage1 (for Channing Cluster): ( Rscript ./test.R ./LV135preX.csv ./LV133postX.csv ) >& testR.log  ##")
message("## Usage2 (for Channing Cluster): qbR ./test.R ./LV135preX.csv ./LV133postX.csv                       ##")
message("## Result: Usage 1 works; Usage 2 (qbR command line) failed to pass the arguments                    ##")
message("## > args[1]                                                                                        ##")
message("## [1] NA                                                                                          ##")
message("## > args[2]                                                                                      ##")
message("## [1] NA                                                                                        ##")
message("## > testfile1 <- args[1]                                                                       ##")
message("## > testfile2 <- args[2]                                                                      ##")
message("## > testfile1                                                                                ##")
message("## [1] NA                                                                                    ##")
message("## > testfile2                                                                              ##")
message("## [1] NA                                                                                  ##")
message("## > LV135preX <- read.csv(testfile1, stringsAsFactors=FALSE)                             ##")
message("## Error in file(file, 'rt') : cannot open the connection                                ##")
message("##########################################################################################")

R.version
R.home()
.libPaths()
sessionInfo()
# installed.packages()
# .libPaths("/udd/retwc/R/library/3.1")
# .libPaths( c( .libPaths(), "/udd/retwc/R/library/3.1") )
# detach("package:GGtools" , unload=TRUE)
# install.packages("name-of-your-package", lib="~/R/library")

message("###############################################################################")
message("# Environmental Path Variables+Working directory settings, R packages loading #")
message("###############################################################################")
# require(reshape2, lib.loc="/udd/retwc/R/library/3.1")
# require("peer", lib.loc="/udd/retwc/R/library/3.1")
# require(reshape2)
# require(MatrixEQTL)
# require(peer)
# setwd("/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522") # for Channing Cluster
# setwd("/n/data1/hms/genetics/seidman/danny/RNASEQ/alextwc/20150630PEER") # for HMS Orchestra Cluster
setwd("/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/20150730PEER")
# setwd("D:/BWHMS/WeeklyMeetings/2015-07-13")
# bestK <- 10 # according to eQTL test using SNDA data
require("VariantAnnotation")
require("GenomicRanges")
require("peer", lib.loc="/udd/retwc/R")
require("GGtools")
source("./InputFiles/pifdrNA.R") # load and execute a script of R commands
vcfpath_01 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh1.vcf-IR.vcf.gz"
vcfpath_02 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh2.vcf-IR.vcf.gz"
vcfpath_03 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh3.vcf-IR.vcf.gz"
vcfpath_04 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh4.vcf-IR.vcf.gz"
vcfpath_05 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh5.vcf-IR.vcf.gz"
vcfpath_06 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh6.vcf-IR.vcf.gz"
vcfpath_07 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh7.vcf-IR.vcf.gz"
vcfpath_08 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh8.vcf-IR.vcf.gz"
vcfpath_09 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh9.vcf-IR.vcf.gz"
vcfpath_10 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh10.vcf-IR.vcf.gz"
vcfpath_11 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh11.vcf-IR.vcf.gz"
vcfpath_12 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh12.vcf-IR.vcf.gz"
vcfpath_13 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh13.vcf-IR.vcf.gz"
vcfpath_14 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh14.vcf-IR.vcf.gz"
vcfpath_15 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh15.vcf-IR.vcf.gz"
vcfpath_16 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh16.vcf-IR.vcf.gz"
vcfpath_17 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh17.vcf-IR.vcf.gz"
vcfpath_18 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh18.vcf-IR.vcf.gz"
vcfpath_19 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh19.vcf-IR.vcf.gz"
vcfpath_20 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh20.vcf-IR.vcf.gz"
vcfpath_21 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh21.vcf-IR.vcf.gz"
vcfpath_22 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh22.vcf-IR.vcf.gz"
vcfpath_23 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh23.vcf-IR.vcf.gz"
vcfpath_24 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh24.vcf-IR.vcf.gz"
vcfpath_25 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh25.vcf-IR.vcf.gz"
vcfpath_26 <- "/proj/regeps/regep00/studies/TRANSCRIBE/analyses/retwc/eQTL_20150522/VCF-IR/Merged142FinalReportCh26.vcf-IR.vcf.gz"
tf01 = TabixFile(vcfpath_01)
tf02 = TabixFile(vcfpath_02)
tf03 = TabixFile(vcfpath_03)
tf04 = TabixFile(vcfpath_04)
tf05 = TabixFile(vcfpath_05)
tf06 = TabixFile(vcfpath_06)
tf07 = TabixFile(vcfpath_07)
tf08 = TabixFile(vcfpath_08)
tf09 = TabixFile(vcfpath_09)
tf10 = TabixFile(vcfpath_10)
tf11 = TabixFile(vcfpath_11)
tf12 = TabixFile(vcfpath_12)
tf13 = TabixFile(vcfpath_13)
tf14 = TabixFile(vcfpath_14)
tf15 = TabixFile(vcfpath_15)
tf16 = TabixFile(vcfpath_16)
tf17 = TabixFile(vcfpath_17)
tf18 = TabixFile(vcfpath_18)
tf19 = TabixFile(vcfpath_19)
tf20 = TabixFile(vcfpath_20)
tf21 = TabixFile(vcfpath_21)
tf22 = TabixFile(vcfpath_22)
tf23 = TabixFile(vcfpath_23)
tf24 = TabixFile(vcfpath_24)
tf25 = TabixFile(vcfpath_25)
tf26 = TabixFile(vcfpath_26)

args <- commandArgs(TRUE)
exprX <- read.csv(args[1], stringsAsFactors=FALSE)
load(args[2]) # The wide fomat of covariance table that prepared by previous R script
if(file.exists(".RData")) load(".RData") else{
    message("# Please load expression dataset...")}

message("###############################################################################")
message("# Picking the entire 2439 genes of Chr20~Chr23 from expr dataset at this step #")
message("# But it became 1410 genes left after the filter of 'Gene Expression has to   #")
message("# be larger than 0.1RPKM from more than 80% subjects expression data (GTEx)'  #")
message("###############################################################################")
# LV135preX <- read.csv("./LV135preX.csv", header=T, check.names=FALSE, stringsAsFactors=FALSE)
exprX$ugene_id <- sub('X$', '', exprX$ugene_id) # remove the postfixed "X" protection
rownames(exprX) <- exprX[,1]
expr <- exprX[, -1]
rm(exprX)
# LVpreChr20_22 <- LV135pre[rownames(LV135pre) %in% rownames(subset(LV135pre, chr=="chr20" | chr=="chr21" | chr=="chr22")), ]
# LVpreChr20_22 <- subset(LV135pre, chr=="chr20" | chr=="chr21" | chr=="chr22") # 23079 genes --> 1474 genes left
# LVpreSubChr <- subset(LV135pre, chr=="chr1" | chr=="chr20" | chr=="chr21" | chr=="chr22") # 23079 genes --> 3787 genes left
# LVpreSubChr <- subset(LV135pre, chr=="chr2" | chr=="chr20" | chr=="chr21" | chr=="chr22") # 23079 genes --> 2928 genes left
# LVpreSubChr <- subset(LV135pre, chr=="chr19" | chr=="chr20" | chr=="chr21" | chr=="chr22") # 23079 genes --> 3077 genes left
exprChr20_23 <- subset(expr, chr=="chr20" | chr=="chr21" | chr=="chr22" | chr=="chr23") # 23079 genes --> 2439 genes left
exprChr20_23 <- exprChr20_23[,-c(1:4)]
exprChr20_23 <- exprChr20_23[!(rowSums(exprChr20_23==0.0) == ncol(exprChr20_23)),] # 2439 genes --> 2267 genes left
exprChr20_23 <- exprChr20_23[rowSums(exprChr20_23>0.1) >= (ncol(exprChr20_23)*0.8),] # 2267 genes --> 1410 genes left
# colnames(exprChr20_23) <- gsub(".*_(B.*V)", "\\1", colnames(exprChr20_23)) # change BxxxxV_BxxxxV to BxxxxV
# colnames(exprChr20_23) <- gsub("(B.*V)", "\\1_\\1", colnames(exprChr20_23)) # change BxxxxV to BxxxxV_BxxxxV
message("# Is there any Missing Value (NA) in expression dataset ???")
any(is.na(exprChr20_23))
message("# Perform Log2 transformation [log2(1+RPKM)] to normalize the dataset")
exprChr20_23 <- log2(exprChr20_23+1)
message("# Saving the current dataset after low-level processing")
save(exprChr20_23, file="./exprChr20_23.rda")
d <- gsub("-", "_", Sys.Date())
write.csv(exprChr20_23, file=paste("./exprChr20_23_", d, ".csv", sep=""), row.names=TRUE)
rm(expr)

message("###############################################################################")
message("# Preparing Gene Coordinates and phenotypes for regressOut()                  #")
message("###############################################################################")
mergedID <- read.table("./InputFiles/geneid.mergedID", header=TRUE, stringsAsFactors=FALSE)
mergedID$ugene_id <- make.names(mergedID$ugene_id, unique=TRUE)
phenotypes <- read.csv("./InputFiles/Covariances.csv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
# phenotypes <- subset(phenotypes, exclude==0, select=-c(ID,Sequence,Reads,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
phenotypes <- subset(phenotypes, exclude==0, select=-c(AXC_time,ID,Sequence,Reads,LVEF,Weight,Height,BMI,CAD,DM,HgbA1C,BAV,CKMB_day1,beta_blocker,exclude))
phenotypes$Subject_id <- gsub("(B.*V)", "\\1_\\1", phenotypes$Subject_id) # change BxxxxV to BxxxxV_BxxxxV
rownames(phenotypes) <- phenotypes[, 1]
phenotypes <- phenotypes[,-1]
exprChr20_23 <- exprChr20_23[, intersect(rownames(phenotypes), colnames(exprChr20_23))]
phenotypes <- phenotypes[intersect(rownames(phenotypes), colnames(exprChr20_23)), ]
save(phenotypes, file="./phenotypes.rda")
# d <- Sys.Date()
# d <- gsub("-", "_", Sys.Date())
write.csv(phenotypes, file=paste("./phenotypes_", d, ".csv", sep=""), row.names=TRUE)

message("###############################################################################")
message("# Run eQTL with PEER k factors plus the known covariance                      #")
message("# 1. run PEER with covariates, e.g. PEER_setCovariates(model, as.matrix(Covs))#")
message("# 2. get residuals from PEER, e.g. residuals = PEER_getResiduals(model),      #")
message("# 3. use the acquired residuals as new expression input                       #")
message("#    for GGtools with or without regressOut(Age, Sex, SeqCenter)              #")
message("###############################################################################")
message("# Run PEER with k factors & eQTL by GGTools...")
# require(peer)
# load("./covs.rda")

# K=c(0,1)
K=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,55,60,70,80,90,100)
for(k in K){
    message("# Starting the PEER factors analyses with k = ", k)
    model = PEER()
    PEER_setCovariates(model, as.matrix(Covs)) # Include known covariates (NxC matrix) into the model
    PEER_setPhenoMean(model, as.matrix(t(exprChr20_23))) # Include the expression dataset t(GxN matrix) into the model
    # PEER_setNk(model, 0)  # no hidden factors
    # PEER_setAdd_mean(model, TRUE)
    PEER_setNk(model,k) # 135/4=33.75 , so I choose "n_unobserved_factors=35"
    getNk <- PEER_getNk(model)
    message("# Starting the iterations with PEER_getNk = ", getNk)
    PEER_setAdd_mean(model, TRUE)
    PEER_setNmax_iterations(model, 10000000)
    PEER_update(model)
    # X = PEER_getX(model)  # now, the getX() factors should include known covariates + mean
    # W = PEER_getW(model)
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = t(PEER_getResiduals(model))  # convert to GxN
    # residuals = PEER_getResiduals(model) # the code that Patrick Evans used which did not transpose the matrix
    rownames(factors) = rownames(Covs)
    rownames(weights) = rownames(exprChr20_23)
    rownames(residuals) = rownames(exprChr20_23)
    colnames(residuals) = colnames(exprChr20_23)
    pdf(file=paste("./diagnostics_peer_", k, ".exprChr20_23.pdf", sep=""), paper="usr", width = 0, height = 0)
    PEER_plotModel(model)
    mtext(paste0("PEER factor = ", k), side=3, line=1)
    dev.off()
    # write.table(residuals_k, file="./residuals_PEER_k.LV135pre.xls", sep="\t")
    # write.table(factors_k, file="./factors_PEER_k.LV135pre.xls", sep="\t")
    # write.table(weights_k, file="./weights_PEER_k.LV135pre.xls", sep="\t")
    write.csv(residuals, file=paste("./residuals_PEER_", k, ".exprChr20_23.csv", sep=""), row.names=TRUE)
    write.csv(factors,   file=paste("./factors_PEER_", k, ".exprChr20_23.csv", sep=""), row.names=TRUE)
    write.csv(weights,   file=paste("./weights_PEER_", k, ".exprChr20_23.csv", sep=""), row.names=TRUE)
    write.csv(precision,   file=paste("./precision_PEER_", k, ".exprChr20_23.csv", sep=""), row.names=FALSE)
    pdf(file=paste("./precision_peer_", k, ".exprChr20_23.pdf", sep=""))
    plot(precision)
    mtext(paste0("PEER factor = ", k), side=3, line=-1)
    dev.off()
    # expr <- residuals + apply(exprChr20_23, 1, mean) 
    # add rowMean to each gene expression residual value; 1 means applying the function "mean" on each row; 2 means on each column
    expr <- residuals
    rm(model)
    rm(factors)
    rm(precision)
    rm(residuals)
    rm(weights)
    OBJ <- ls()
    message("# Here are the currently existed objects: "); print(OBJ)
    message("# This is the cycle K = ", k)
    message("#   ")
    message("###############################################################################")
    message("# Making SE (SummarizedExperiment) object then performing regressOut()        #")
    message("# by using GenomicRanges package                                              #")
    message("###############################################################################")
    message("# Starting making SE(SummarizedExperiment) objects with or without regressOut()")
    expr <- as.data.frame(expr)
    expr$ugene_id <- rownames(expr)
    names(expr)[names(expr)=="expr$ugene_id"] <- "ugene_id"
    expr <- merge(mergedID, expr, by="ugene_id")
    expr$chr <- gsub("chrX", "chr23", expr$chr)
    expr$chr <- gsub("chrY", "chr24", expr$chr)
    rowD <- GRanges(seqnames=expr$chr, IRanges(expr$start, expr$end), strand=expr$strand)
    assay <- data.matrix(expr[, -c(1:5)])
    exprSE <- SummarizedExperiment(assay, rowRanges=rowD, colData=DataFrame(phenotypes))
    rownames(exprSE) = expr$ugene_id
    colnames(exprSE) = rownames(phenotypes)
    save(exprSE, file=paste("./exprSE_PEER_", k, ".exprChr20_23.rda", sep=""))
    # corrExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter)+Ischemic_time+ReadsLength)
    corrExprSE = regressOut(exprSE, ~Sex+Age+factor(SeqCenter))
    save(corrExprSE, file=paste("./corrExprSE_PEER_", k, ".exprChr20_23.rda", sep=""))
    message("# Finished making SE(SummarizedExperiment) objects and saved them as exprSE(without RegressOut) and corrExprSE(with regressOut)")
    
    message("###############################################################################")
    message("# Run GGTools with Genotyping VCF files on corrExprSE object                  #")
    message("###############################################################################")
    # require(GGtools)
    message("# Starting running GGtools package")
    message("# This is still the cycle K = ", k)
    message("#   ")
    message("# Now is running GGtools with corrExprSE(with regressOut) objects")
    message("# Checking for universal heterozygous loci for exclusion")
    corrExprSE20 = corrExprSE[ seqnames(corrExprSE) == "chr20", ]
    print(table(seqnames(corrExprSE20[1:10,])))
    seqlevelsStyle(corrExprSE20) = "NCBI"
    run20 = cisAssoc( corrExprSE20, tf20, cisradius=100000, lbmaf=.05 )
    chisq = run20$chisq
    perm1 = run20$permScore_1
    perm2 = run20$permScore_2
    perm3 = run20$permScore_3
    fullFDR_run20 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run20 <- run20[which(!is.na(run20$chisq))] 
    print(length(fullFDR_run20))
    print(length(run20))
    run20$fdr = fullFDR_run20
    save(run20, file=paste("./run20.corrExprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    corrExprSE21 = corrExprSE[ seqnames(corrExprSE) == "chr21", ]
    print(table(seqnames(corrExprSE21[1:10,])))
    seqlevelsStyle(corrExprSE21) = "NCBI"
    run21 = cisAssoc( corrExprSE21, tf21, cisradius=100000, lbmaf=.05 )
    chisq = run21$chisq
    perm1 = run21$permScore_1
    perm2 = run21$permScore_2
    perm3 = run21$permScore_3
    fullFDR_run21 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run21 <- run21[which(!is.na(run21$chisq))] 
    print(length(fullFDR_run21))
    print(length(run21))
    run21$fdr = fullFDR_run21
    save(run21, file=paste("./run21.corrExprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    corrExprSE22 = corrExprSE[ seqnames(corrExprSE) == "chr22", ]
    print(table(seqnames(corrExprSE22[1:10,])))
    seqlevelsStyle(corrExprSE22) = "NCBI"
    run22 = cisAssoc( corrExprSE22, tf22, cisradius=100000, lbmaf=.05 )
    chisq = run22$chisq
    perm1 = run22$permScore_1
    perm2 = run22$permScore_2
    perm3 = run22$permScore_3
    fullFDR_run22 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run22 <- run22[which(!is.na(run22$chisq))] 
    print(length(fullFDR_run22))
    print(length(run22))
    run22$fdr = fullFDR_run22
    save(run22, file=paste("./run22.corrExprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    corrExprSE23 = corrExprSE[ seqnames(corrExprSE) == "chr23", ]
    print(table(seqnames(corrExprSE23[1:10,])))
    seqlevelsStyle(corrExprSE23) = "NCBI"
    run23 = cisAssoc( corrExprSE23, tf23, cisradius=100000, lbmaf=.05 )
    chisq = run23$chisq
    perm1 = run23$permScore_1
    perm2 = run23$permScore_2
    perm3 = run23$permScore_3
    fullFDR_run23 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run23 <- run23[which(!is.na(run23$chisq))] 
    print(length(fullFDR_run23))
    print(length(run23))
    run23$fdr = fullFDR_run23
    save(run23, file=paste("./run23.corrExprSE_PEER_", k, ".rda", sep=""))
    run20_fdr.05 <- run20[which(run20$fdr<=0.05)]
    run21_fdr.05 <- run21[which(run21$fdr<=0.05)]
    run22_fdr.05 <- run22[which(run22$fdr<=0.05)]
    run23_fdr.05 <- run23[which(run23$fdr<=0.05)]
    df.run20_fdr.05 <- as(run20_fdr.05, "data.frame")
    df.run21_fdr.05 <- as(run21_fdr.05, "data.frame")
    df.run22_fdr.05 <- as(run22_fdr.05, "data.frame")
    df.run23_fdr.05 <- as(run23_fdr.05, "data.frame")
    corrExprSeQTL_fdr.05 <- rbind(df.run20_fdr.05, df.run21_fdr.05, df.run22_fdr.05, df.run23_fdr.05)
    save(corrExprSeQTL_fdr.05, file=paste("./corrExprSeQTL_fdr.05_PEER_", k, ".rda", sep=""))
    print(length(corrExprSeQTL_fdr.05))
    print(dim(corrExprSeQTL_fdr.05))
    assPairs <- dim(corrExprSeQTL_fdr.05)
    eGenes <- nlevels(factor(corrExprSeQTL_fdr.05$probeid))
    message("# The current PEER factor = ", k)
    message("# The number of SNP-Gene association pairs have FDR <= 5% in regressOut() group is ", assPairs[1])
    message("# The count of cis-eQTL genes with regressOut() is ", eGenes)
    corrExprSeQTL_fdr.05$eGenesCount <- nlevels(factor(corrExprSeQTL_fdr.05$probeid))
    write.csv(corrExprSeQTL_fdr.05, file=paste("./corrExprSeQTL_fdr.05_", k, ".csv", sep=""), row.names=TRUE)
    
    message("#   ")
    message("# Now is running GGtools with exprSE(without regressOut) objects")
    message("# Checking for universal heterozygous loci for exclusion")
    exprSE20 = exprSE[ seqnames(exprSE) == "chr20", ]
    print(table(seqnames(exprSE20[1:10,])))
    seqlevelsStyle(exprSE20) = "NCBI"
    run20 = cisAssoc( exprSE20, tf20, cisradius=100000, lbmaf=.05 )
    chisq = run20$chisq
    perm1 = run20$permScore_1
    perm2 = run20$permScore_2
    perm3 = run20$permScore_3
    fullFDR_run20 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run20 <- run20[which(!is.na(run20$chisq))] 
    print(length(fullFDR_run20))
    print(length(run20))
    run20$fdr = fullFDR_run20
    save(run20, file=paste("./run20.exprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    exprSE21 = exprSE[ seqnames(exprSE) == "chr21", ]
    print(table(seqnames(exprSE21[1:10,])))
    seqlevelsStyle(exprSE21) = "NCBI"
    run21 = cisAssoc( exprSE21, tf21, cisradius=100000, lbmaf=.05 )
    chisq = run21$chisq
    perm1 = run21$permScore_1
    perm2 = run21$permScore_2
    perm3 = run21$permScore_3
    fullFDR_run21 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run21 <- run21[which(!is.na(run21$chisq))] 
    print(length(fullFDR_run21))
    print(length(run21))
    run21$fdr = fullFDR_run21
    save(run21, file=paste("./run21.exprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    exprSE22 = exprSE[ seqnames(exprSE) == "chr22", ]
    print(table(seqnames(exprSE22[1:10,])))
    seqlevelsStyle(exprSE22) = "NCBI"
    run22 = cisAssoc( exprSE22, tf22, cisradius=100000, lbmaf=.05 )
    chisq = run22$chisq
    perm1 = run22$permScore_1
    perm2 = run22$permScore_2
    perm3 = run22$permScore_3
    fullFDR_run22 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run22 <- run22[which(!is.na(run22$chisq))] 
    print(length(fullFDR_run22))
    print(length(run22))
    run22$fdr = fullFDR_run22
    save(run22, file=paste("./run22.exprSE_PEER_", k, ".rda", sep=""))
    message("# Checking for universal heterozygous loci for exclusion")
    exprSE23 = exprSE[ seqnames(exprSE) == "chr23", ]
    print(table(seqnames(exprSE23[1:10,])))
    seqlevelsStyle(exprSE23) = "NCBI"
    run23 = cisAssoc( exprSE23, tf23, cisradius=100000, lbmaf=.05 )
    chisq = run23$chisq
    perm1 = run23$permScore_1
    perm2 = run23$permScore_2
    perm3 = run23$permScore_3
    fullFDR_run23 = pifdrNA(chisq, c(perm1,perm2,perm3))
    run23 <- run23[which(!is.na(run23$chisq))] 
    print(length(fullFDR_run23))
    print(length(run23))
    run23$fdr = fullFDR_run23
    save(run23, file=paste("./run23.exprSE_PEER_", k, ".rda", sep=""))
    run20_fdr.05 <- run20[which(run20$fdr<=0.05)]
    run21_fdr.05 <- run21[which(run21$fdr<=0.05)]
    run22_fdr.05 <- run22[which(run22$fdr<=0.05)]
    run23_fdr.05 <- run23[which(run23$fdr<=0.05)]
    df.run20_fdr.05 <- as(run20_fdr.05, "data.frame")
    df.run21_fdr.05 <- as(run21_fdr.05, "data.frame")
    df.run22_fdr.05 <- as(run22_fdr.05, "data.frame")
    df.run23_fdr.05 <- as(run23_fdr.05, "data.frame")
    exprSeQTL_fdr.05 <- rbind(df.run20_fdr.05, df.run21_fdr.05, df.run22_fdr.05, df.run23_fdr.05)
    save(exprSeQTL_fdr.05, file=paste("./exprSeQTL_fdr.05_PEER_", k, ".rda", sep=""))
    print(length(exprSeQTL_fdr.05))
    print(dim(exprSeQTL_fdr.05))
    assPairs <- dim(exprSeQTL_fdr.05)
    eGenes <- nlevels(factor(exprSeQTL_fdr.05$probeid))
    message("# The current PEER factor = ", k)
    message("# The number of SNP-Gene association pairs have FDR <= 5% in without-regressOut() group is ", assPairs[1])
    message("# The count of cis-eQTL genes without regressOut() is ", eGenes)
    exprSeQTL_fdr.05$eGenesCount <- nlevels(factor(exprSeQTL_fdr.05$probeid))
    write.csv(exprSeQTL_fdr.05, file=paste("./exprSeQTL_fdr.05_", k, ".csv", sep=""), row.names=TRUE)
    message("#  ")
    message("###############################################################################")
}

timing <- proc.time()
message("# Finishing The entire R Script for time = ")
print(timing)

