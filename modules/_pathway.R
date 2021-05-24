#!/usr/bin/env Rscript
###########################################
# R script for running pathway analysis
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 2/21/2021
# version: 1.0
# Usage: Rscript $0 -i DEresult.all.CONDITION_treated_vs_control.xls.gz -F medium --genome rn6 -o pathway/FDR05
# Input: input file should include columns: 
# - symbol (for GO/GSEA/SPIA)
# - padj (for medium and strigent filter), 
# - log2FoldChange (for SPIA, GSEA)
# - lfcSE (for strigent filter)
# - stat (for GSEA)
# Note: Upgrate to use msigdb_v7.2 (2/20/2021)
###########################################
# install packages
require("optparse",quietly=T, warn.conflicts=F) || install.packages('optparse', repo='http://cran.revolutionanalytics.com');

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path for input file from DE analysis result", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="Output directory path [default=%default]", metavar="character"),
  make_option(c("-N", "--topN"), type="integer", default=20, 
              help="Top N result selected [default=%default]", metavar="integer"),
  make_option(c("-q", "--qcutoff"), type="double", default=0.05, 
              help="Adjust p-value cutoff [default=%default]", metavar="double"),
  make_option(c("-F", "--filters"), type="character", default='medium', 
              help="Filter of input data [default=%default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default='hg19', 
              help="Genome assembly [default=%default]", metavar="character"),
  make_option(c("-p", "--program"), type="character", default='all', 
              help="GO|GSEA|SPIA [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Input result file from DE analysis must be supplied", call.=FALSE)
}

suppressPackageStartupMessages(library('tidyverse',logical.return=T)) || {install.packages('tidyverse'); suppressPackageStartupMessages(library('tidyverse',logical.return=T))}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
suppressPackageStartupMessages(library('SPIA',logical.return=T)) || {BiocManager::install('SPIA',ask=F); suppressPackageStartupMessages(library('SPIA',logical.return=T))}
suppressPackageStartupMessages(library('org.Hs.eg.db',logical.return=T)) || {BiocManager::install("org.Hs.eg.db",ask=F); suppressPackageStartupMessages(library('org.Hs.eg.db',logical.return=T))}
suppressPackageStartupMessages(library('biomaRt',logical.return=T)) || {BiocManager::install('biomaRt',ask=F); suppressPackageStartupMessages(library('biomaRt',logical.return=T))}
suppressPackageStartupMessages(library("topGO",logical.return=T)) || {BiocManager::install("topGO",ask=F); suppressPackageStartupMessages(library("topGO",logical.return=T))}
suppressPackageStartupMessages(library('fgsea',logical.return=T)) || {BiocManager::install('fgsea',ask=F); suppressPackageStartupMessages(library('fgsea',logical.return=T))}
suppressPackageStartupMessages(library('pathview',logical.return=T)) || {BiocManager::install('pathview',ask=F); suppressPackageStartupMessages(library('pathview',logical.return=T))}

input_DE_result=opt$input
output_dir=opt$out
topN=opt$topN
Q_CUTOFF=opt$qcutoff
index=opt$genome
filters=opt$filters
program=opt$program

if(is.null(program)) program="all"

message(paste("-i", input_DE_result,"-o",output_dir,"-N",topN,"-q",Q_CUTOFF,"-g",index,"-F",filters))

# debug:
# input_DE_result="DEresult.all.CONDITION_treated_vs_control.xls.gz"; index='rn6'; topN=20;Q_CUTOFF=0.05; filters='medium'; output_dir="~/neurogen/rnaseq_Rot/results/merged/DE"
# input_DE_result="DEresult.padj_05.Cell_type_DP_vs_DN_in_HC.xls"; index='hg19'; topN=20;Q_CUTOFF=0.05; filters='medium'; output_dir="~/rnaseq_MS/results/DE"
# input_DE_result="Merge_circexp_norm_filtered_and_enriched.2genes.stat.xls"; index='hg19'; topN=20;Q_CUTOFF=0.05;filters='no';  output_dir="."
# input_DE_result="DEresult.DE_SNDA.CONDITION2_PD_vs_HC.xls"; index='hg19'; topN=20; Q_CUTOFF=0.05;filters='no';  output_dir="~/projects/circRNA/results/DE_SNDA"
# input_DE_result="DEresult.DE.rn6.RotClen_vs_RotVeh.all.CONDITION_RotClen_vs_RotVeh.xls.gz"; index='rn6'; topN=20;Q_CUTOFF=0.05; filters='medium'; output_dir="~/neurogen/rnaseq_Rot/results/merged/DE.rn6.RotClen_vs_RotVeh"
# setwd("~/projects/circRNA/results/DE2gene_SNDA"); input_DE_result="DEresult.DE2gene_SNDA.CONDITION2_PD_vs_HC.xls"; index='hg19'; topN=20; Q_CUTOFF=0.05;filters='loose';  output_dir="pathway"
# setwd("~/projects/andrew2020/results/DE"); input_DE_result="DEresult.padj05_log2FCgt1.CONDITION_media_vs_mock.xls"; index='hg19'; topN=20;Q_CUTOFF=0.05; filters='medium'; output_dir="pathway"； program="all"


# check input
if(!file.exists(input_DE_result)) {stop(paste(input_DE_result, "doesn't exist. Exit!"), call.=FALSE);}

file.exists(output_dir) || dir.create(output_dir, recursive = T)

genome_name=switch(index, "hg19" = "Homo_sapiens", "mm9" = "Mus_musculus", "rn6" = "Rattus_norvegicus")
gname=tolower(sub("(.).*_(.*)","\\1\\2",genome_name))  # e.g. Homo_sapiens --> hsapiens
GENOME_DIR=file.path("~/neurogen/referenceGenome", genome_name, "UCSC", index)

###################################################
message("#Step1: loading data")
###################################################
if(grepl("DiagnosisMS_vs_HC",input_DE_result)) DE_result=read.table(input_DE_result,check.names =F, stringsAsFactors=F) else 
  DE_result=read.delim(input_DE_result, stringsAsFactors=F, row.names = 1, header=T, check.names =F)

if(nrow(DE_result)<=5) {stop(paste(input_DE_result, "has too few (<=5) records. Exit!"), call.=FALSE);}
if(!("symbol" %in% colnames(DE_result))) {
  message("The input file does not have a 'symbol' colume.")
  if("geneName" %in% colnames(DE_result)) {
    DE_result$symbol = DE_result$geneName
    message("Use 'geneName' colume instead.")
  } else{
    message("Column names are:")
    message(paste(colnames(DE_result), collapse = "\n"))
    cat("Enter the column name for gene symbol:")
    symbol_col_name <- readLines("stdin",n=1)
    if(!(symbol_col_name %in% colnames(DE_result))) stop("The input file must include 'symbol' in the colume. Exit!", call.=FALSE);
    DE_result$symbol = DE_result[[symbol_col_name]]
  }
}

# get human ortholog if the DE result is from non-human species
DE_result$human_orth_symbol=DE_result$symbol
if(index!='hg19') {
  #ensembl = useEnsembl(biomart="ensembl", version=90); grep("rnorvegicus",listDatasets(ensembl)$dataset);message("done");quit('no');
  gmart = useEnsembl("ensembl", version=90, dataset = paste0(gname, "_gene_ensembl")) # version 90 is the first version with rat_gene_ensembl in biomart (9/17/2018)
  #listAttributes(gmart, page='homologs'); listFilters(gmart);
  attributes = c("ensembl_gene_id",
                 "hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_perc_id", "hsapiens_homolog_perc_id_r1")
  orth.human = getBM(attributes,
                     filters="with_hsapiens_homolog",values=TRUE, mart = gmart,
                     uniqueRows=TRUE)
  #dim(orth.human);
  orth.human = mutate(orth.human, ave_perc_id=hsapiens_homolog_perc_id+hsapiens_homolog_perc_id_r1) %>% 
    group_by(ensembl_gene_id) %>% 
    filter(ave_perc_id==max(ave_perc_id)) %>% 
    filter(1:n() == 1) # make sure each gene only assigned one ortholog (the best recipocal one)
  DE_result$human_orth_symbol=orth.human$hsapiens_homolog_associated_gene_name[match(rownames(DE_result), orth.human$ensembl_gene_id)] # the unfound will be NA
}
#head(DE_result)
dim(DE_result); DE_result=DE_result[!is.na(DE_result$human_orth_symbol),]; dim(DE_result)

# to save all genes for GSEA
DE_result0=DE_result;

# filter
if(filters == 'no') DE_result=DE_result # all i.e. p<=1
if(filters == 'loose') DE_result=subset(DE_result, pvalue<=.05) # pvalue<0.05
if(filters == 'medium') DE_result=subset(DE_result, padj<=.05)
if(filters == 'strigent') DE_result=subset(DE_result, padj<=0.05 & lfcSE<2)

if(nrow(DE_result)<=5) {stop(paste(input_DE_result, "has too few (<=5) records. Exit!"), call.=FALSE);}

setwd(output_dir)


if(("symbol" %in% colnames(DE_result)) && (program=="all" || grepl("GO",program))) {
  ####################################################
  message("#Step2: running GO terms enrichment anlaysis")
  ####################################################
  ## require: only gene symbol
  
  source("~/neurogen/pipeline/RNAseq/bin/lib.R")
  df <- read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F)
  gene_names = df$V7[df$V8=="protein_coding"]
  
  # all DE genes
  topGOenrichment(DE_result$human_orth_symbol, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all', output=file.path(input_DE_result))
  if("log2FoldChange" %in% colnames(DE_result)){
    # UP regulated DE genes
    topGOenrichment(DE_result$human_orth_symbol[DE_result$log2FoldChange>0], allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='up', output=file.path(input_DE_result))
    # DOWN regulated DE genes
    topGOenrichment(DE_result$human_orth_symbol[DE_result$log2FoldChange<0], allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='down', output=file.path(input_DE_result))
  }
  
  ## simply over-representation analysis (ORA)
  ORA(inputGenes=DE_result$human_orth_symbol, allGenes=gene_names, nCutoff=3, pCutoff=0.01, output=file.path(input_DE_result))
  if("log2FoldChange" %in% colnames(DE_result)){
    # UP regulated DE genes
    ORA(DE_result$human_orth_symbol[DE_result$log2FoldChange>0], allGenes=gene_names, output=file.path(paste0(input_DE_result,".up")))
    # DOWN regulated DE genes
    ORA(DE_result$human_orth_symbol[DE_result$log2FoldChange<0], allGenes=gene_names, output=file.path(paste0(input_DE_result,".down")))
  }
}

if(all(c('symbol', 'stat') %in% colnames(DE_result0)) && (program=="all" || grepl("GSEA",program))){
  ###################################################
  message("#Step3: GSEA pathway analysis")
  ###################################################
  ## require: gene symbol and stat
  ## require no cutoff for the input, e.g. use all genes 
  
  # Extract the symbol and stat (as metric) 
  ranks = DE_result0 %>% dplyr::select(symbol=human_orth_symbol, stat=stat) %>%
    group_by(symbol) %>% top_n(n=1, wt=abs(stat)) %>% ungroup()
  
  # OR use -log10(p-value)*sign(2logFC) as metric  (08312018)
  # see discussion on the metric: https://www.biostars.org/p/159029/
  ranks <- data.frame(symbol =  DE_result0$human_orth_symbol, stat = -log10(DE_result0$pval) * sign(DE_result0$log2FoldChange))
  
  ranks = setNames(ranks$stat, ranks$symbol)
  
  gt=data.frame();
  pdf(file.path(paste0(input_DE_result, '.fGSEA.pdf')))
  for(i in c("c2.cp","c2.cp.kegg", "c2.cgp","c3.tft", "c7.all")){
    message(i);
    
    gmt.file = paste0("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/msigdb_v7.2/msigdb_v7.2_GMTs/",i,".v7.2.symbols.gmt")
    
    pathways = gmtPathways(gmt.file)
    fgseaRes = fgsea(pathways, stats=ranks, minSize=10, maxSize=500, nPermSimple=10000)
    fgseaRes = fgseaRes[order(fgseaRes$padj),]  # default is BH-adjusted p-value, same as FDR
    
    # get FDR < 0.05 pathways
    fgseaRes = subset(fgseaRes, padj < Q_CUTOFF)

    if(nrow(fgseaRes)>0) {
      # plot
      for(x in fgseaRes$pathway) {p=plotEnrichment(pathways[[x]], stats=ranks) + labs(title=x); print(p);}
      
      # flatten the list of leadingEdge
      fgseaRes$leadingEdge = vapply(fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
      
      gt=rbind(gt, data.frame(gene_set=i, fgseaRes))
    }
  }
  dev.off()
  write.table(gt, file=file.path(paste0(input_DE_result, '.fGSEA.xls')), sep="\t", quote=F, row.names=F, col.names=T)
  
} 

## require: only gene symbol and log2FoldChange
if(all(c('symbol', 'log2FoldChange') %in% colnames(DE_result)) && (program=="all" || grepl("SPIA",program))){
  ###################################################
  message("#Step4: SPIA analysis")
  ###################################################
  ## Prerequisition
  # download all KEGG pathway and install in SPIA
  # cd /Library/Frameworks/R.framework/Versions/4.0/Resources/library/SPIA
  # curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,":"); print "curl http://rest.kegg.jp/get/"a[2]"/kgml -o extdata/keggxml/hsa/"a[2]".xml"}' | bash
  # curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,":"); print "curl http://rest.kegg.jp/get/"a[2]"/image -o extdata/keggxml/hsa/"a[2]".png"}' | bash
  #library(SPIA)
  #makeSPIAdata(kgml.path=system.file("extdata/keggxml/hsa",package="SPIA"),organism="hsa",out.path="./extdata")
  
  
  # add ENTREZ ID using the human_orth_symbol (as Pathway data might be more enriched in human data)
  # Note: ortholog genes (e.g. PAX6 in human and Pax6 in mouse) have different ENTREZ ID
  human_orth_Entrez = mapIds(org.Hs.eg.db, 
                             keys = DE_result$human_orth_symbol, 
                             column = "ENTREZID", 
                             keytype = "SYMBOL", 
                             multiVals='first')
  DE_result$human_orth_Entrez = as.numeric(human_orth_Entrez) # NA for those genes without Entrez ID
  degGIs <- as.vector(DE_result$log2FoldChange)
  names(degGIs) <- DE_result$human_orth_Entrez
  # remove NA
  degGIs <- degGIs[!is.na(names(degGIs))]
  # remove duplicate
  degGIs=degGIs[!duplicated(names(degGIs))]
  
  res <- spia(de=degGIs, 
              all=mappedkeys(org.Hs.egGENENAME), 
              organism='hsa', 
              nB=2000, plots=FALSE, 
              beta=NULL, 
              combine="fisher")
  
  ## A data frame containing the ranked pathways and various statistics: 
  # pSize is the number of genes on the pathway; 
  # NDE is the number of DE genes per pathway; 
  # tA is the observed total preturbation accumulation in the pathway; 
  # pNDE is the probability to observe at least NDE genes on the pathway using a hypergeometric model; 
  # pPERT is the probability to observe a total accumulation more extreme than tA only by chance; 
  # pG is the p-value obtained by combining pNDE and pPERT; 
  # pGFdr and pGFWER are the False Discovery Rate and respectively Bonferroni adjusted global p-values; 
  # Status gives the direction in which the pathway is perturbed (activated or inhibited). 
  # KEGGLINK gives a web link to the KEGG website that displays the pathway image with the differentially expressed genes highlighted in red.
  
  # add color code to the KEGGLINK based on fold-change of each DE gene 
  # Potential issue: some pathways are not in the Pathview web server and it will return an error (see https://support.bioconductor.org/p/90506/)
  library(pathview) # BiocManager::install("pathview",ask=F);
  sapply(filter(res, pGFdr<Q_CUTOFF, !(ID %in% c("04723","04320","04215","05206")), file.exists(paste0("/Library/Frameworks/R.framework/Versions/4.0/Resources/library/SPIA/extdata/keggxml/hsa/hsa",ID,".xml")), !file.exists(paste0("hsa",ID,".",input_DE_result,".png"))) %>% dplyr::select(ID), 
         function(pid) pathview(gene.data  = degGIs, 
                                pathway.id = pid, 
                                species = "hsa", 
                                kegg.dir = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/SPIA/extdata/keggxml/hsa/",
                                out.suffix = input_DE_result,
                                limit  = list(gene=max(abs(degGIs)), cpd=1)))
  
  # make text file and plot
  write.table(subset(res, pGFdr<Q_CUTOFF), file.path(paste0(input_DE_result, '.SPIA.xls')), sep="\t", quote=F, row.names=F, col.names=T)
  
  pdf(file.path(paste0(input_DE_result, '.SPIA.pdf')))
  plotP(res, threshold=Q_CUTOFF)
  res$ID=res$Name;
  plotP(res, threshold=Q_CUTOFF)
  dev.off()
} 