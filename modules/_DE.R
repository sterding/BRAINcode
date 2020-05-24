###########################################
# A general framework in R for running differntial expression analysis using DESeq2
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 5/23/2018
# version: 1.0
# Usage: Rscript --vanilla ../../src/_DE.R -i genes.htseqcount.cufflinks.mm9.uniq.xls -c covariate.mm9.txt -O m,i -g mm9 -o DE.mm9
# Requirement:
# 1. The headers in covariate file should be in CAPITAL format for SAMPLE_ID (required), SUBJECT_ID (if any), CELL_TYPE (optional), CONDITION (required)
# TODO: The current code is only for comparison in levels of CELL_TYPE and/or CONDITION. In future, we want to compare any levels the user assigned.
###########################################
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input expression dataset file name", metavar="character"),
  make_option(c("-c", "--covariate"), type="character", default=NULL, 
              help="Covariate file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output directory path [default=%default]", metavar="character"),
  make_option(c("-O", "--output_addition_columns"), type="character", default=NULL, 
              help="output additioal columns than DEseq2 default", metavar="character"),
  make_option(c("-G", "--genome"), type="character", default='hg19', 
              help="Genome assembly [default=%default]", metavar="character"),
  make_option(c("-r", "--reference_level"), type="character", default='control', 
              help="Reference level for comparison [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_expression_filename=opt$input
input_covariance_filename=opt$covariate
output_dir=opt$out
output_additonal_columns=opt$output_addition_columns
index=opt$genome
ref=opt$reference_level

if(is.null(input_expression_filename) || is.null(input_covariance_filename)){
  print_help(opt_parser)
  stop("Both input expression matrix and covariate file must be supplied", call.=FALSE)
}

# debug
# setwd("~/neurogen/rnaseq_PD/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.TCPY.uniq.xls"; input_covariance_filename="TCPY.covariates.txt"; output_dir="TCPY"; index="hg19"
# setwd("~/neurogen/rnaseq_Rot/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.mm9.uniq.xls"; input_covariance_filename="covariate.mm9.txt"; output_dir="DE.mm9"; index="mm9"; output_additonal_columns="mi"
# setwd("~/neurogen/rnaseq_Rot/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.rn6.uniq.xls"; input_covariance_filename="covariate.rn6.txt"; output_dir="DE.rn6"; index="rn6"; output_additonal_columns="mi"
# setwd("~/neurogen/rnaseq_Rot/results/merged"); input_expression_filename="genes.htseqcount.cufflinks.rn6.uniq.xls"; input_covariance_filename="covariate.rn6.txt2"; output_dir="DE.rn6.VehClen_vs_VehVeh"; index="rn6"; output_additonal_columns="mi"

# check input
if(!file.exists(input_expression_filename)) {stop(paste(input_expression_filename, "doesn't exist. Exit!"), call.=FALSE);}
if(!file.exists(input_covariance_filename)) {stop(paste(input_covariance_filename, "doesn't exist. Exit!"), call.=FALSE);}
genome_name=switch(index, "hg19" = "Homo_sapiens", "mm9" = "Mus_musculus", "rn6" = "Rattus_norvegicus")

# Create folder if the directory doesn't exist
dir.create(file.path(output_dir,'report/figures'), recursive =T, showWarnings = FALSE)

pwd=getwd()
setwd(output_dir)

# install packages
require('tidyverse',quietly=T, warn.conflicts=F) || install.packages('tidyverse', repo='http://cran.revolutionanalytics.com');
require('RCurl',quietly=T, warn.conflicts=F) || install.packages('RCurl', repo='http://cran.revolutionanalytics.com');
require('hexbin',quietly=T, warn.conflicts=F) || install.packages('hexbin', repo='http://cran.revolutionanalytics.com');
require('pheatmap',quietly=T, warn.conflicts=F) || install.packages('pheatmap', repo='http://cran.revolutionanalytics.com');
require('RColorBrewer',quietly=T, warn.conflicts=F) || install.packages('RColorBrewer', repo='http://cran.revolutionanalytics.com');
require('hwriter',quietly=T, warn.conflicts=F) || install.packages('hwriter', repo='http://cran.revolutionanalytics.com');
source("https://bioconductor.org/biocLite.R"); 
require('vsn',quietly=T, warn.conflicts=F) || biocLite('vsn');
require('DESeq2',quietly=T, warn.conflicts=F) || biocLite('DESeq2');
#require('ReportingTools',quietly=T, warn.conflicts=F) || biocLite('ReportingTools');
require('BiocParallel',quietly=T, warn.conflicts=F) || biocLite('BiocParallel');
require('limma',quietly=T, warn.conflicts=F) || biocLite('limma');

comparisons_colors=c("green","red","orange","purple","blue","pink")

###########################################
message("#step1: load data...")
###########################################
if(file.exists("DESeq2.RData")) load("DESeq2.RData") else {
  
  # debug
  # genome_name="Rattus_norvegicus"; index='rn6'; input_expression_filename='genes.htseqcount.cufflinks.allSamples.uniq.xls'; input_covariance_filename='covariate.txt'
  #annotation
  #genes_annotation = read.table(file.path("~/neurogen/referenceGenome",genome_name,"UCSC",index,"Annotation/Genes/gencode.annotation.genes.bed+2"), header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType"));
  genes_annotation = read.table(file.path("~/neurogen/referenceGenome",genome_name,"UCSC",index,"Annotation/Genes/annotation.genes.bed6+3"), sep="\t", quote="", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType","geneDescription"));
  # remove the tailing version number in the ENSEMBL ID, if any
  genes_annotation$geneID = sub("\\..*","",genes_annotation$geneID)
  
  # raw reads count
  if(tolower(tools::file_ext(input_expression_filename) == "rds")) {
    cts=readRDS(file.path(pwd,input_expression_filename))
  }else cts=read.delim(file.path(pwd,input_expression_filename), row.names = 1,check.names =F)
  # remove those non-geneID rows, e.g. __no_feature (pre-mRNA reads) and __ambiguous (see http://htseq.readthedocs.io/en/master/count.html )
  dim(cts); cts=cts[grep("^__", rownames(cts), invert = T),]; dim(cts);
  
  # covariance table
  covarianceTable = read.table(file.path(pwd,input_covariance_filename), sep="\t", header = T,check.names =F)
  head(covarianceTable)
  rownames(covarianceTable) = covarianceTable$SAMPLE_ID;
  
  # subset and re-order
  all(rownames(covarianceTable) %in% colnames(cts))
  dim(cts); cts = cts[, intersect(colnames(cts), rownames(covarianceTable))]; dim(cts)
  dim(covarianceTable); covarianceTable = covarianceTable[colnames(cts), ]; dim(covarianceTable);
  covarianceTable[] <- lapply(covarianceTable, function(x) if(is.factor(x)) factor(x) else x) # drop levels  after subsetting
  all(rownames(covarianceTable) == colnames(cts))
  
  # check the structure of corariance table, to see if it's necessary to factorize some of the columns
  str(covarianceTable)
  
  # # factorize
  # covarianceTable = mutate(covarianceTable,
  #                          CONDITION=factor(CONDITION, levels = c("control","treated")))
  # str(covarianceTable)
  
  # relevel
  covarianceTable$CONDITION = relevel(covarianceTable$CONDITION, ref=ref)
  
  levels(covarianceTable$CONDITION)
  
  
  ###########################################
  message("#step2: load data to DEseq")
  ###########################################
  
  # Note: With no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level.
  # Ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  
  # making formula
  fmla <- names(covarianceTable)[!names(covarianceTable) %in% c("SAMPLE_ID", "SUBJECT_ID", "Age","PMI")]
  fmla <- as.formula(paste(" ~ ", paste(fmla, collapse= "+")))
  fmla
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = covarianceTable,
                                design= fmla)
  
  ## pre-filtering
  #keep <- rowSums(counts(dds)) >= 10
  keep <- rowSums(counts(dds) >= 5) >= 4  # at least 4 samples with more than 5 reads
  dim(dds); dds <- dds[keep,]; dim(dds);
  #head(sort(rowMeans(counts(dds)), decreasing = T))
  # head(counts(dds)[order(rowMeans(counts(dds)), decreasing =T),])
  
  save.image(file.path("DESeq2.RData"))

}
# script to generate html report page
makeNewImages <- function(df,...){
  imagename <- c()
  tablename <- c()
  for (i in 1:nrow(df)){
    ensId <- rownames(df)[i]
    symbol <- df$symbol[i]    
    imagename[i] <- paste('plot', ensId, symbol, 'pdf', sep = ".")
    tablename[i] <- paste('plot', ensId, symbol, 'txt', sep = ".")
    
    d <- data.frame(samples=colnames(assay(ntd)), 
                    expression_ntd=assay(ntd)[ensId,], 
                    expression_raw=assay(dds)[ensId,], 
                    condition=colData(dds)$CONDITION)
    N=1
    if("CELL_TYPE" %in% names(colData(dds))) {
      d=cbind(d,cell=colData(dds)$CELL_TYPE)
      N=length(levels(colData(dds)$CELL_TYPE))
      p=ggplot(d, aes(x=condition, y=expression_ntd, group=cell)) + 
        geom_boxplot(position=position_dodge(.8), width=.8, outlier.shape = NA) +
        geom_jitter(size=1.5, position = position_jitter(width=.15)) +
        facet_wrap(~ cell, nrow=1) +
        stat_summary(fun.y=mean, geom="line", colour="red", size=0.8) +
        xlab("CONDITION") + ylab("log2(counts+1)") + ggtitle(paste(ensId, "(",symbol,")"))
    } else {
      N=length(levels(colData(dds)$CONDITION))
      p=ggplot(d, aes(x=condition, y=expression_ntd)) + 
        geom_boxplot(position=position_dodge(.8), width=.5, outlier.shape = NA) +
        geom_jitter(size=1.5, position = position_jitter(width=.15)) +
        theme_bw() +
        xlab("CONDITION") + ylab("log2(counts+1)") + ggtitle(symbol,subtitle =ensId)
    }
    if(!file.exists(file.path('report/figures',tablename[i]))) {
      write.table(d,file.path('report/figures',tablename[i]),sep="\t", quote =F, row.names=F, col.names = T)
    }
    if(!file.exists(file.path('report/figures',imagename[i]))) {
      #png(file.path('report/figures', imagename[i]), height = 250, width = 600)
      pdf(file.path('report/figures', imagename[i]), height = 4, width = 3*N)
      print(p)
      dev.off()
    }
  }
  ## Using the following code to show thumb figures. It's slow to display if many
  # df$Boxplot <- hwriteImage(paste0('figures/', imagename), 
  #                           link=paste0('figures/', imagename), 
  #                           table=FALSE, width=100)
  df$Boxplot <- hwrite('boxplot', link = paste0('figures/', imagename), table=F)
  df$Rawdata <- hwrite("data", link = paste0('figures/', tablename), table=F)
  df$symbol <- hwrite(as.character(df$symbol), 
                      link = paste0("http://useast.ensembl.org/",genome_name,"/Gene/Summary?db=core;g=",as.character(rownames(df))), 
                      table=F)
  return(df)
}
###########################################
message("#step3: QA of the data [optional]")
###########################################
# Note: This part is not necessary for DEseq, but important for data QA

##--------------------------------------
## 3.1: compare different vairance stablization methods
##--------------------------------------

ntd <- normTransform(dds) # log2(x+1)
vsd <- vst(dds, blind=T) # Note: blind to the design, equal to design = ~ 1

# # using limma to remove covariates, it returns adjusted values in log2 scale
# vsd_adjusted_log2 <- removeBatchEffect(assay(vsd), batch=vsd$Batch, batch2=vsd$Sex, covariates = colData(vsd)[,c('Age','PMI')])

pdf("diagnosis.pdf")
par(mar=c(10,5,3,3));boxplot(assay(ntd), las=2, cex.axis=.7, ylab="log2(x+1) transform")
msd <- meanSdPlot(counts(dds), ranks = FALSE); msd$gg + ggtitle("no transformation")
msd <- meanSdPlot(assay(ntd), ranks = FALSE); msd$gg + ggtitle("log2(x+1) transform")
msd <- meanSdPlot(assay(vsd), ranks = FALSE); msd$gg + ggtitle("VST")
#msd <- meanSdPlot(vsd_adjusted_log2, ranks = FALSE); msd$gg + ggtitle("vsd_adjusted_log2")
dev.off()

##--------------------------------------
## 3.2: save normalized reads count
##--------------------------------------

## save the raw reads count 
write.table(counts(dds), 
            file=paste0(input_expression_filename, ".filtered.raw.xls"), 
            sep="\t", quote = F, 
            col.names = NA, row.names = TRUE)
## save the variance-stabilized data
write.table(assay(vsd), 
            file=paste0(input_expression_filename, ".filtered.vsd.xls"), 
            sep="\t", quote = F, 
            col.names = NA, row.names = TRUE)
# ## save the vsd adjusted
# write.table(vsd_adjusted_log2, 
#             file=paste0(input_expression_filename, ".filtered.vsd_adjusted_log2.xls"), 
#             sep="\t", quote = F, 
#             col.names = NA, row.names = TRUE)

##--------------------------------------
## 3.3: clustering of samples
##--------------------------------------

pdf("clustering.tree.pdf")

#sampleDists <- as.dist((1 - cor(assay(vsd))))
sampleDists <- dist(t(assay(vsd)))  # scale on columns; dist on rows

## heatmap
par(cex=0.5, mar=c(5, 8, 4, 1))

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         main = paste0(output_dir,": heatmap"),
         fontsize = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         scale = "row", 
         clustering_method = 'ward.D', 
         cluster_rows = TRUE,
         clustering_distance_rows = 'euclidean',
         cluster_cols = T,
         clustering_distance_cols = 'euclidean')

## tree
plot(hclust(sampleDists,method = "ward.D"), xlab='', main=paste0(output_dir,": Cluster Dendrogram"))

## PCA
pcaData <- plotPCA(vsd, intgroup = c( "CONDITION"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p=ggplot(pcaData, aes(x = PC1, y = PC2, color = CONDITION)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(paste0(output_dir,": PCA")) + 
  coord_fixed()  
print(p)

dev.off()

###########################################
message("#step4: Run DE")
###########################################

##################################################################################################
##################################################################################################
# If neither CONDITION or CELL_TYPE is present, exit with a warning message
##################################################################################################
##################################################################################################
if(sum(c("CONDITION", "CELL_TYPE") %in% names(covarianceTable))==0){
  stop("No CONDITION or CELL_TYPE detected. Please check your covariate file!", call. =F);
}
##################################################################################################
##################################################################################################
# If CONDITION or CELL_TYPE are present, not both, call DE individually (no subsetting, no interaction)
##################################################################################################
##################################################################################################
if(sum(c("CONDITION", "CELL_TYPE") %in% names(covarianceTable))==1){
  fmla <- names(covarianceTable)[!names(covarianceTable) %in% c("SAMPLE_ID", "SUBJECT_ID", "Age","PMI")]
  design(dds) <- as.formula(paste(" ~ ", paste(fmla, collapse= "+")))
  
  # In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
  dds <- DESeq(dds, betaPrior=T, parallel=TRUE, BPPARAM=MulticoreParam(4))  
  resultsNames(dds)
  
  comparisons=data.frame(stringsAsFactors=F);
  if("CONDITION" %in% names(covarianceTable)) comparisons=rbind(comparisons, cbind("CONDITION", t(combn(levels(covarianceTable$CONDITION),2))))
  if("CELL_TYPE" %in% names(covarianceTable)) comparisons=rbind(comparisons, cbind("CELL_TYPE", t(combn(levels(covarianceTable$CELL_TYPE),2))))
  
  apply(comparisons, 1, function(x){
    com_name= paste(x[1], x[3], "vs", x[2], sep="_")
    message(paste("processing comparison:",com_name))
    res <- results(dds, contrast = c(x[1], x[3], x[2]),   # contract format: factor name, numerator in the fold change, denominator in the fold change
                   alpha = 0.05, 
                   parallel=TRUE, BPPARAM=MulticoreParam(4))
    ## Shrink log2 fold changes [required for version >=1.16]
    #res2 <- lfcShrink(dds, contrast = c(x[1], x[2], x[3]), res=res, type='normal', parallel=TRUE, BPPARAM=MulticoreParam(4))
    ## You can get the shrunken LFC either with lfcShrink like above or with betaPrior=TRUE. It will be the same shrunken LFC and the same as previously calculated in DESeq2. The difference is that betaPrior=TRUE will give you a p-value for the shrunken LFC, while lfcShrink (at the moment) is only giving you the LFC, and is keeping the p-value for the test of the MLE LFC. 
    ## see https://support.bioconductor.org/p/95695/ and https://support.bioconductor.org/p/98833/#98843
    
    summary(res)
    head(res); dim(res)
    # decimal value of Fold-change
    res$FoldChange <- 2**res$log2FoldChange
    # add annotation
    res$symbol <- genes_annotation$geneSymbol[match(sub("\\..*","",row.names(res)), genes_annotation$geneID)]
    res$geneType <- genes_annotation$geneType[match(sub("\\..*","",row.names(res)), genes_annotation$geneID)]
    res$geneDescription <- genes_annotation$geneDescription[match(sub("\\..*","",row.names(res)), genes_annotation$geneID)]
    
    # add additional columns in the output
    if(!is.null(output_additonal_columns) && grepl("m", output_additonal_columns)){ # mean of raw expression values for each group, 
      if("CONDITION" %in% names(covarianceTable)) baseMeanPerLvl <- sapply( levels(dds$CONDITION), function(lvl) rowMeans( counts(dds,normalized=FALSE)[,dds$CONDITION == lvl] ) )
      if("CELL_TYPE" %in% names(covarianceTable)) baseMeanPerLvl <- sapply( levels(dds$CELL_TYPE), function(lvl) rowMeans( counts(dds,normalized=FALSE)[,dds$CELL_TYPE == lvl] ) )
      colnames(baseMeanPerLvl) = paste0("baseMean_raw.", colnames(baseMeanPerLvl))
      res = cbind(res, baseMeanPerLvl)
    }
    if(!is.null(output_additonal_columns) && grepl("M", output_additonal_columns)){ # mean of normalized expression values for each group, 
      if("CONDITION" %in% names(covarianceTable)) baseMeanPerLvl <- sapply( levels(dds$CONDITION), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$CONDITION == lvl] ) )
      if("CELL_TYPE" %in% names(covarianceTable)) baseMeanPerLvl <- sapply( levels(dds$CELL_TYPE), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$CELL_TYPE == lvl] ) )
      colnames(baseMeanPerLvl) = paste0("baseMean_norm.", colnames(baseMeanPerLvl))
      res = cbind(res, baseMeanPerLvl)
    }
    if(!is.null(output_additonal_columns) && grepl("i", output_additonal_columns)){ # individual raw expression values of each sample 
      individual <- counts(dds,normalized=FALSE)
      colnames(individual) = paste0("ind_raw.", colnames(individual))
      res = cbind(res, individual)
    }
    if(!is.null(output_additonal_columns) && grepl("I", output_additonal_columns)){ # individual normalized expression values of each sample
      individual <- counts(dds,normalized=TRUE)
      colnames(individual) = paste0("ind_norm.", colnames(individual))
      res = cbind(res, individual)
    }
    res <- res[order(res$padj),]
    head(res); dim(res)
    
    # write to xls
    # write.table(as.data.frame(subset(res, padj<0.05)), 
    #             file=file.path(paste0("DEresult.padj_05.", com_name ,".xls")), 
    #             sep="\t", quote =F, na="", row.names=T, col.names = NA)
    # write.table(as.data.frame(subset(res, pvalue<0.05)), 
    #             file=file.path(paste0("DEresult.pvalue_05.", com_name ,".xls")), 
    #             sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    write.table(as.data.frame(res), 
                file=file.path(paste0("DEresult.",output_dir,".all.", com_name ,".xls.gz")), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    write.table(as.data.frame(subset(res, padj<=0.05 & abs(log2FoldChange)>=1)), 
                file=file.path(paste0("DEresult.",output_dir,".padj05_log2FCgt1.", com_name ,".xls")), 
                sep="\t", quote =F, na="", row.names=T, col.names = NA)
    
    # ## Note: 10% of non-significant(NS) genes are randomly selected in order to increase the performance of the generating graph
    NS=subset(res, padj>0.05 | abs(log2FoldChange)<1)
    n_NS=nrow(NS)
    NS=NS[sample(n_NS,round(n_NS*.20)),]
    dim(res); res=DESeqResults(rbind(NS, subset(res, padj<=0.05 & abs(log2FoldChange)>=1))); dim(res);
    
    ## MAKING PLOTS
    pdf(file.path(paste0("DEresult.",output_dir,".padj_05.", com_name ,".pdf")), paper = 'USr')
    # scatter plot
    # TOADD: 
    
    ## ==============================
    # MA plot
    ## ==============================
    DESeq2::plotMA(res, alpha = 0.05, colNonSig = "gray", main=paste0(output_dir,": ",com_name))
    
    ## ==============================
    # vocano plot
    ## ==============================
    topT <- as.data.frame(res)
    with(topT, plot(log2FoldChange, -log10(padj), 
                    pch=20, cex=0.5, main=paste0(output_dir,": ",com_name), col='gray',
                    xlab=bquote(~Log[2]~fold~change), 
                    ylab=bquote(~-log[10]~FDR)))
    if(nrow(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1))>0){
      with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1), 
         points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1))
      with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1), 
          text(log2FoldChange, -log10(padj), labels=symbol, cex=0.5, pos=1, offset=0.2))
    }
    abline(v=c(0,-1, 1), lty=c(3,4,4), lwd=c(1,2,2))
    abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
    
    
    if(nrow(subset(res, abs(log2FoldChange)>=1 & padj<0.05))>0){    
      ## ==============================
      message("# heatmap for top 10 DE genes")
      ## ==============================
      topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
      
      #topT=topT[order(-topT$baseMean),] # sort by baseMean in decreasing order
      #topT=topT[order(-abs(topT$log2FoldChange)),] # sort by abs(log2FoldChange) in decreasing order
      topT=topT[order(topT$padj),] # sort by padj in increasing order
      topT=rbind(head(subset(topT, log2FoldChange<0),10),head(subset(topT, log2FoldChange>0),10)) # top 10 down-regulated and top 10 up-regulated
      topDE=assay(vsd)[rownames(topT),]
      rownames(topDE) = topT$symbol
      annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
      rownames(annotation_row) = topT$symbol
      annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
      ann_colors = list(
        updown = c(up = "red", down = "blue"),
        geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue'),
        CONDITION = c("green","red")
      )
      names(ann_colors$CONDITION)=c(x[3],x[2])
      
      
      ## Scale/center each genes (by rows)
      topDE=t(scale(t(as.matrix(topDE))))
      ## trim max and min
      MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
      
      par(cex=0.5, mar=c(5, 8, 4, 1))
      pheatmap(topDE,
               fontsize = 8,
               main =paste0(output_dir,": heatmap for top 10 DE genes"),
               #fontsize_col = 5,
               #width=8, height=5,filename=paste0("DEresult.padj_05.", com_name ,".heatmap.pdf"),
               border_color = NA,
               color = colorRampPalette(c("blue", "white", "red"))(50),
               annotation_row = annotation_row,
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               drop_levels = TRUE,
               scale = "none", 
               clustering_method = 'ward.D', 
               cluster_rows = TRUE,
               clustering_distance_rows = "correlation",
               cutree_rows = 2,cutree_cols=2,
               cluster_cols = TRUE,
               clustering_distance_cols = "correlation")
      
      ## ==============================
      message("# heatmap for top 20 DE genes")
      ## ==============================
      topT=subset(res, abs(log2FoldChange)>=1 & padj<0.05)
      topT=topT[order(topT$padj),] # sort by padj in increasing order
      topT=rbind(head(subset(topT, log2FoldChange<0),20),head(subset(topT, log2FoldChange>0),20)) # top 10 down-regulated and top 10 up-regulated
      topDE=assay(vsd)[rownames(topT),]; rownames(topDE) = topT$symbol
      annotation_row = dplyr::select(as.data.frame(topT), geneType, log2FoldChange, symbol) %>% mutate(updown=ifelse(log2FoldChange>0,"up","down")) %>% select(updown, geneType)
      rownames(annotation_row) = topT$symbol
      annotation_col = dplyr::select(as.data.frame(colData(dds)), CONDITION)
      ann_colors = list(
        updown = c(up = "red", down = "blue"),
        geneType = c(protein_coding = "darkblue", lincRNA = "orange", Mt_rRNA='pink', antisense='yellow', pseudogene='gray', processed_transcript='lightgray',Mt_tRNA='purple',misc_RNA='red',snoRNA='lightblue'),
        CONDITION = c("green","red")
      )
      names(ann_colors$CONDITION)=c(x[3],x[2])
      
      ## Scale/center each genes (by rows)
      topDE=t(scale(t(as.matrix(topDE))))
      ## trim max and min
      MAX=2; MIN=-2; topDE[topDE>MAX]=MAX;topDE[topDE<MIN]=MIN
      
      par(cex=0.5, mar=c(5, 8, 4, 1))
      pheatmap(topDE,
               main =paste0(output_dir,": heatmap for top 20 DE genes"),
               fontsize = 8,
               #fontsize_col = 5,
               border_color = NA,
               color = colorRampPalette(c("blue", "white", "red"))(50),
               annotation_row = annotation_row,
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               drop_levels = TRUE,
               scale = "none", 
               clustering_method = 'ward.D', 
               cluster_rows = TRUE,
               clustering_distance_rows = "correlation",
               cutree_rows = 2,cutree_cols=2,
               cluster_cols = TRUE,
               clustering_distance_cols = "correlation")
      
      
    #   message("## Exporting results to HTML")
    #   htmlRep <- HTMLReport(shortName=com_name, title=com_name,
    #                         reportDirectory="report")
    #   publish(makeNewImages(head(as.data.frame(subset(res, abs(log2FoldChange)>=1 & padj<0.05, select=!grepl("^ind_", colnames(res)))),1000)), htmlRep)
    #   finish(htmlRep)
    # 
      }
    
    dev.off()     
  })
}

## scp html result to web server
system(paste0("rsync -a --chmod=u+rwx,g+rwx,o+rwx report/* xd010@panda.dipr.partners.org:~/public_html/DE_reports/",output_dir))
message(paste0("Visit the HTML result at: http://panda.partners.org/~xd010/DE_reports/",output_dir))
message("Done!")