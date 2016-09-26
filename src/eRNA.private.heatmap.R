##############################################
## R script to create heatmpa for the private eRNAs from three major cell types
## Author: Xianjun Dong
## Usage: Rscript $HOME/neurogen/pipeline/RNAseq/src/eRNA.private.heatmap.R ~/eRNAseq/eRNA.private.major.merged.bed ~/eRNAseq/eRNA.private.major.merged.meanRPM.allSamples.xls 
## Date:
## Version: 0.0
############################################

args<-commandArgs(TRUE)
geneAnnotation=args[1]  
geneExpression=args[2]

geneAnnotation='~/eRNAseq/eRNA.private.major.merged.bed'
geneExpression='~/eRNAseq/eRNA.private.major.merged.meanRPM.allSamples.xls'

setwd("~/eRNAseq")

# ================================
## read data
# ================================

covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples
require(RCurl)
covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
covs=subset(covs, BRAINCODE.final.selection==1)

geneAnnotation=read.delim(geneAnnotation, header=F, stringsAsFactors=F)
colnames(geneAnnotation)=c('chr','start','end','ID','private_HTNE')

geneExpression=read.delim(geneExpression, header = T, check.names =F)
#head(geneExpression)
rownames(geneExpression)=geneExpression[,1]; 
geneExpression=geneExpression[,covs$sampleName];
exp=geneExpression*1000 + .1  # meanRPM -> RPKM
#colnames(exp)=gsub(".*_.*_(.*)_.*_rep.*", "\\1", colnames(exp))
plot(hclust(as.dist(1-cor(exp)), method="ward.D"), cex=0.5)

# reorder the columns
colorder=c("SNDA","TCPY","MCPY","PBMC", "FB")
exp=exp[,unlist(lapply(colorder, grep, colnames(exp)))]
# median by group
exp_group=data.frame(SNDA=apply(exp[,grep("SNDA",colnames(exp))], 1, median),
                PY=apply(exp[,grep("PY",colnames(exp))], 1, median),
                NN=apply(exp[,grep("FB|PBMC",colnames(exp))], 1, median),
                private_HTNE=geneAnnotation$private_HTNE[match(rownames(exp), geneAnnotation$ID)])
exp_group$private_HTNE=factor(exp_group$private_HTNE, levels = c("SNDAonly","PYonly","NNonly"))
head(exp_group)

# calcluate fold_change: for each cell type X, fold_change = X / max(not X)
exp_group = cbind(exp_group, fc=apply(exp_group, 1, function(x) {
    ct=sub("only","",x[4]); 
    y=as.numeric(x[1:3]);
    if(ct=="SNDA") fc=y[1]/max(y[2],y[3]);
    if(ct=="PY") fc=y[2]/max(y[1],y[3]);
    if(ct=="NN") fc=y[3]/max(y[2],y[1]);
    fc
}))

# Option1: take the top 100 private TNEs with highest median expression in each group 
# sort by private_HTNE, then by expression of each group
sn=subset(exp_group, private_HTNE=='SNDAonly'); 
py=subset(exp_group, private_HTNE=='PYonly'); 
nn=subset(exp_group, private_HTNE=='NNonly'); 
topN=100
exp_group_top = rbind(head(sn[order(-sn$SNDA),],topN), head(py[order(-py$PY),],topN), head(nn[order(-nn$NN),],topN))

# Option2: take the private genes with fc>=30
# sort by fc in each group
#exp_group = exp_group[with(exp_group, order(private_HTNE, -fc)),]
#exp_group_top=subset(exp_group, fc>=30) # N=2978
#table(exp_group_top$private_HTNE)

# # Option3: take the top 10% of private genes with the most FC
# library(dplyr)
# exp_group_top = cbind(name=rownames(exp_group), exp_group) %>%   # add rownames as a column as dplyr::filter discard rownames
#     group_by(private_HTNE) %>%
#     arrange(private_HTNE, desc(fc)) %>%
#     filter(fc > quantile(fc, .95))
# exp_group_top=as.data.frame(exp_group_top)
# rownames(exp_group_top)=exp_group_top[,1]; exp_group_top=exp_group_top[,-1]

dim(exp_group_top)
table(exp_group_top$private_HTNE)
tapply(exp_group_top$fc, exp_group_top$private_HTNE, summary)
# $SNDAonly
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 457.3   574.7   718.1   867.4   993.9  5286.0 
# 
# $PYonly
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1360    1546    1828    2062    2288   11620 
# 
# $NNonly
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 192.0   620.0   857.1   972.4  1179.0  5195.0 

exp_selected = exp[rownames(exp_group_top),] 

#exp_selected=10^exp_selected

# scale each row to [0,1]
exp_selected_scaled = (exp_selected - apply(exp_selected, 1, min)) / (apply(exp_selected, 1, max) - apply(exp_selected, 1, min))
# trim a bit
exp_selected_scaled[exp_selected_scaled>0.85]=1
#exp_selected_scaled[exp_selected_scaled<0.15]=0

# scale and center each row to u(0,1)
# exp_selected_scaled = t(scale(t(exp_selected))) # scale is for the columns
# exp_selected_scaled[exp_selected_scaled>2]=2
# exp_selected_scaled[exp_selected_scaled<-2]=-2

range(exp_selected_scaled)

plot(hclust(as.dist(1-cor(exp_selected_scaled)), method="ward.D2"), cex=0.5)

# visualize it
# source("http://bioconductor.org/biocLite.R"); biocLite(c("RColorBrewer", "pheatmap"))
library("pheatmap")
library("RColorBrewer")

# Generate annotations
cnames=sub("(.*rep[0-9]).*","\\1",colnames(exp_selected_scaled))
colnames(exp_selected_scaled)=cnames
annotation_col=data.frame(do.call(rbind, strsplit(cnames,"_")), stringsAsFactors=F)
rownames(annotation_col) = cnames
colnames(annotation_col) = c("Diagnosis","SubjectID","Cell_Type","Batch")
annotation_col=annotation_col[,"Cell_Type", drop=F]
Cell_Type2=annotation_col$Cell_Type
Cell_Type2[grep("FB|PBMC",annotation_col$Cell_Type)]="NN"
Cell_Type2[grep("PY",annotation_col$Cell_Type)]="PY"
Cell_Type2=factor(Cell_Type2, levels = c('SNDA','PY','NN'))
annotation_col$Cell_Type=factor(annotation_col$Cell_Type, levels = colorder)

annotation_row=geneAnnotation[match(rownames(exp_selected_scaled), geneAnnotation$ID),'private_HTNE', drop=F]
rownames(annotation_row)= rownames(exp_selected_scaled)
annotation_row$private_HTNE=factor(annotation_row$private_HTNE, levels = c('SNDAonly','PYonly','NNonly'))

cCell_Type = c(SNDA="#F22A7B", TCPY="#3182bd", MCPY='#2659B2', PBMC='#513931', FB='#BC9371'); 
cPrivacy = c(SNDAonly="#F22A7B", PYonly="#3182bd", NNonly='#513931'); 
ann_colors = list(Cell_Type = cCell_Type, private_HTNE = cPrivacy)

bk=seq(min(exp_selected_scaled),max(exp_selected_scaled), length=100)
hmcols<- colorRampPalette(c("black","yellow","yellow"))(length(bk)-1)  

hm.parameters <- list(as.matrix(exp_selected_scaled),
                      breaks=bk,
                      color = hmcols,
                      border_color=NA,
                      fontsize=5,
                      scale = "none",
                      show_rownames = F, show_colnames = F,
                      main = "Heatmap of private HTNEs",
                      cluster_rows = F, 
                      cluster_cols = F,
                      annotation_row = annotation_row,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors,
                      gaps_row=cumsum(table(annotation_row$private_HTNE)),
                      gaps_col=cumsum(table(Cell_Type2)),
                      cellwidth=2,
                      cellheight=.02
)

png(file=paste("private.eRNAs.cluster","png", sep="."), width=600, height=400, res=400)
par(mar=c(0,0,0,0))
rotate <- function(x) t(apply(x, 2, rev)) # function to rotate a matrix with 90 degree clock-wise (ref: http://stackoverflow.com/questions/16496210/rotate-a-matrix-in-r)
image(rotate(as.matrix(exp_selected_scaled)), col=hmcols, breaks=bk, useRaster=F, axes =F)
dev.off()

# To draw to file 
do.call("pheatmap", c(hm.parameters, filename=paste("private.eRNAs.cluster","pdf", sep=".")))
#do.call("pheatmap", c(hm.parameters, filename=paste("private.eRNAs.cluster","png", sep=".")))
## To draw the heat map on screen 
#do.call("pheatmap", hm.parameters)
