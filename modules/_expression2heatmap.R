library(pheatmap) # install.packages('pheatmap')
library(tidyverse)
library(RColorBrewer)
library(optparse)
options(stringsAsFactors=F)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input expression dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output directory path [default=%default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default='hg19', 
              help="Genome assembly [default=%default]", metavar="character"),
  make_option(c("-l", "--gene_list"), type="character", default=NULL, 
              help="The input of list [default=%default]", metavar="character"),
  make_option(c("-t", "--gene_type"), type="character", default='circRNA', 
              help="RNA type [gene|eRNA|circRNA] [default=%default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default='~/projects/circRNA/data/Merge_circexplorer_BC197.filtered.enriched.annotation.bed14.rds', 
              help="circRNA annotation [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_expression_filename=opt$input
output_dir=opt$out
index=opt$genome
annotation_path=opt$annotation
file_of_gene_list = opt$gene_list
gene_type = opt$gene_type

if(is.null(input_expression_filename) || is.null(file_of_gene_list)){
  print_help(opt_parser)
  stop("Both input expression matrix and the list of genes must be supplied", call.=FALSE)
}
if(!file.exists(input_expression_filename)) {stop(paste(input_expression_filename, "doesn't exist. Exit!"), call.=FALSE);}
if(!file.exists(file_of_gene_list)) {stop(paste(file_of_gene_list, "doesn't exist. Exit!"), call.=FALSE);}

# setwd("~/neurogen/rnaseq_PD/results/DE_SNDA.gene"); input_expression_filename="genes.htseqcount.cufflinks.allSamples.BCv2.uniq.xls.filtered.vsd_adjusted_log2.xls"; file_of_gene_list="DEresult.DE_SNDA.gene.MUSS.padj0.05.xls";
# setwd("~/neurogen/rnaseq_Rot/results/merged/DE.rn6.RotClen_vs_RotVeh"); input_expression_filename="genes.htseqcount.cufflinks.rn6.uniq.xls.filtered.vsd.xls"; file_of_gene_list="DEresult.DE.rn6.RotClen_vs_RotVeh.all.CONDITION_RotClen_vs_RotVeh.Mito.padj0.05.xls";

# raw expression file
if(tolower(tools::file_ext(input_expression_filename)) == "rds") {
  exp=readRDS(file.path(input_expression_filename))
}else exp=read.delim(file.path(input_expression_filename), row.names = 1,check.names =F)
head(exp); dim(exp)
## ==============================
message("# heatmap for the input list of genes")
## ==============================
geneList = read.delim(file = file_of_gene_list, row.names = 1, check.names = F, stringsAsFactors = F, header = T)
if(grepl("rnaseq_Rot", getwd())) geneList = rownames_to_column(geneList) %>% mutate(geneID=rowname,geneName=symbol) %>% column_to_rownames()
head(geneList); dim(geneList)

df=exp[rownames(geneList),]; 
rownames(df) = paste(geneList$geneID, " - ", geneList$geneName)

x=cor(t(df), method='spearman')
dim(x)

# # optional: remove row/col which max(abs)<0.7
# x2=x; x2[lower.tri(x2, T)]=0
# INDEX=apply(abs(x2),1,max, na.rm=T)>0.7
# x=x[INDEX,INDEX]

annotation_row = dplyr::select(geneList, one_of(c('geneType', 'log2FoldChange', 'geneName'))) %>% rownames_to_column() %>% 
  rowwise() %>% mutate(updown=ifelse(log2FoldChange>0,"up","down"), rowname=paste(rowname," - ", geneName)) %>% 
  select(rowname, one_of('updown','geneType')) %>% column_to_rownames()

suppressPackageStartupMessages(library('RColorBrewer',logical.return=T) || install.packages('RColorBrewer', repo='http://cran.revolutionanalytics.com'))
suppressPackageStartupMessages(library('pheatmap',logical.return=T) || install.packages('pheatmap', repo='http://cran.revolutionanalytics.com'))

ann_colors = list(
  updown = c(up = "red", down = "blue"))
## add color for geneType
genetypes = sort(table(as.character(geneList$geneType)), decreasing = T)
genetypes = setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(genetypes)), names(genetypes))
ann_colors = c(geneType = list(genetypes), ann_colors)

par(cex=0.5, mar=c(5, 8, 4, 1))
pheatmap(x,
         main =paste0("Co-expression cluster for the gene list:",file_of_gene_list,"\n", "clustering_distance = euclidean, clustering_method = complete"),
         fontsize = 5,
         fontsize_row = 4,
         border_color = NA,
         breaks=seq(-1,1,0.01),
         color = colorRampPalette(c("darkblue", "blue", "white", "red",'darkred'))(200),
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         drop_levels = TRUE,
         cexRow=0.4, labCol ="",
         #clustering_method = 'ward.D',
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         filename = paste0(file_of_gene_list,".heatmap.pdf"),width = nrow(x)/15, height =nrow(x)/15,
         scale = "none",show_colnames =F)
