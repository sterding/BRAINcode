#!/usr/bin/env Rscript
# ===========================================================
# R script to make hivePlot for gene-enhancer correlation
# author: Xianjun Dong
# date: 2016-02-15
# input: matrix of gene expression, enhaner expression, bed file of gene, enhancer position
# Usage: Rscript ~/pipeline/src/correlation2hivePlot.R eRNA.meanRPM.xls /data/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls eRNA.bed /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed chr17_43583680_44506585 
# ===========================================================

args<-commandArgs(TRUE)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("grid")) {
  install.packages("grid", dependencies = TRUE)
  library(grid)
}
if (!require("HiveR")) {
  install.packages("HiveR", dependencies = TRUE)
  library(HiveR)
}

valA=args[1]  # e.g. eRNA.expression.tab
valB=args[2]  # e.g. genes.expression.tab
posA=args[3]  # e.g. eRNA postion bed
posB=args[4]  # e.g. genes postition bed
regions=args[5]
CUTOFF=args[6]

# debug
valA="eRNA.meanRPM.xls"; valB="/data/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls"; 
posA="eRNA.bed"; posB="/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed";
regions="chr4:90,445,404-90,959,294";
CUTOFF=0.4

message("reading data ...")
# ===============================================
valA=read.table(valA, header=T, stringsAsFactors =F)  # e.g. eRNA.expression.tab
rownames(valA) = valA[,1]; valA=valA[,-1]; 
valB=read.table(valB, header=T, stringsAsFactors =F)  # e.g. genes.expression.tab
rownames(valB) = valB[,1]; valB=valB[,-1]; 
colnames(valB)=gsub("FPKM.","",colnames(valB))
colnames(valB)=gsub("_0$","",colnames(valB))

# same order of columns for A and B
common = intersect(colnames(valA), colnames(valB))
valA=valA[,common]; valB=valB[,common]

posA=read.table(posA, header=F, stringsAsFactors =F)  # e.g. eRNA postion bed
colnames(posA) = c('chr','start','end','name') #,'score','strand')
posA$score=0; posA$strand=".";
posB=read.table(posB, header=F, stringsAsFactors =F)  # e.g. genes postition bed
colnames(posB) = c('chr','start','end','name','score','strand','symbol','type')
rownames(posA) = posA$name; rownames(posB) = posB$name

region = strsplit(gsub(",","",regions), "[_:-]")[[1]]
c=region[1]; s=as.numeric(region[2]); e=as.numeric(region[3])

message("calculating the gene-enhancer correlation matrix ...")
# ===============================================
nmsA = subset(posA, chr==c & start<e & end>s, select=name)[,1]
nmsB = subset(posB, chr==c & start<e & end>s, select=name)[,1]
# only show protein-coding genes
#nmsB = subset(posB, chr==c & start<e & end>s & type=='protein_coding', select=name)[,1]
# only MAPT and KANSL1
#nmsB = c('ENSG00000186868.11','ENSG00000120071.8')  
df = cor(t(valA[nmsA,]), t(valB[nmsB,]), method='spearman') # row is A and col is B

# ENSID --> gene symbol
nmsB2symbol = posB[nmsB,'symbol']

# filter abs(cor())<0.5
df[is.na(df)]=0
df[abs(df)<CUTOFF]=0

message("plot the correlation matrix ...")
# ===============================================
# gene-gene interaction
hive1 <- adj2HPD(df, axis.cols =c('#66666655','#66666655'), desc = "Gene-enhancer correaltion")  # now hive1$nodes are first rows and then cols
# set axis
hive1$nodes$axis=as.integer(c(rep(1,nrow(df)),rep(2,ncol(df))))
# set radius of nodes
hive1$nodes$radius=c(ifelse(posA[nmsA, 'strand']=="-",posA[nmsA,'end'],posA[nmsA,'start']), ifelse(posB[nmsB, 'strand']=="-",posB[nmsB,'end'],posB[nmsB,'start']))
minR=min(hive1$nodes$radius, s); maxR=max(hive1$nodes$radius, e)
hive1$nodes$radius=hive1$nodes$radius-minR+100
# set size of nodes
hive1$nodes$size=.4
# set color of nodes
hive1$nodes$color=c(rep('magenta',nrow(df)), rep('black',ncol(df)))
# set color of edge
hive1$edges$color=ifelse(hive1$edges$weight>0, rgb(1,0,0,abs(hive1$edges$weight)), rgb(1,1,1,abs(hive1$edges$weight)))

## add two artificial nodes on each axis, min-100 and max+100
ranger=c(100,maxR-minR+100) # range(hive1$nodes$radius);
hive1$nodes = rbind(hive1$nodes, list(as.integer(max(hive1$nodes$id)+1), '',as.integer(1),ranger[1]-100,0.5,'NA'))
hive1$nodes = rbind(hive1$nodes, list(as.integer(max(hive1$nodes$id)+1), '',as.integer(1),ranger[2]+100,0.5,'NA'))
hive1$nodes = rbind(hive1$nodes, list(as.integer(max(hive1$nodes$id)+1), '',as.integer(2),ranger[1]-100,0.5,'NA'))
hive1$nodes = rbind(hive1$nodes, list(as.integer(max(hive1$nodes$id)+1), '',as.integer(2),ranger[2]+100,0.5,'NA'))
chkHPD(hive1)

# optional: add enhaner-enhancer correlation
x=cor(t(valA[nmsA,]), method='spearman'); x[is.na(x)]=0; x[abs(x)<CUTOFF]=0
h=adj2HPD(x)
h$edges$color=ifelse(h$edges$weight>0, rgb(0,1,0,abs(h$edges$weight)), rgb(1,1,1,abs(h$edges$weight)))
tmp=h$edges$id1; h$edges$id1=h$edges$id2; h$edges$id2=tmp; # reverse the end of eges to flip the curve
hive1$edges = rbind(hive1$edges, h$edges)

# optional: add gene-gene correlation
x=cor(t(valB[nmsB,]), method='spearman'); x[is.na(x)]=0; x[abs(x)<CUTOFF]=0
h=adj2HPD(x)
h$edges$id1=h$edges$id1+nrow(df); h$edges$id2=h$edges$id2+nrow(df); # assign to axis2
#tmp=h$edges$id1; h$edges$id1=h$edges$id2; h$edges$id2=tmp; # reverse the end of eges to flip the curve
h$edges$color=ifelse(h$edges$weight>0, rgb(0,0,1,abs(h$edges$weight)), rgb(1,1,1,abs(h$edges$weight)))
hive1$edges = rbind(hive1$edges, h$edges)

# set weight of edge
hive1$edges$weight=abs(hive1$edges$weight) # or fixed width .5

# remove gene self correlation (e.g. R=1)
hive1 = mineHPD(hive1, option = "remove zero edge")

# normalize edges$weight to [0-1]
#hive1$edges$weight=2*(hive1$edges$weight- min(hive1$edges$weight)) / (max(hive1$edges$weight) - min(hive1$edges$weight)) 

# check
chkHPD(hive1)

# generate nodes.csv file for nodes annotation
# see example by cat(readLines("/PHShome/xd010/R/x86_64-unknown-linux-gnu-library/3.1/HiveR/extdata/Misc/HECticks.txt"), sep="\n")
nodes=subset(hive1$nodes, axis==2 & lab!="", select=lab)
nodes=cbind(nodes, text=posB[nodes$lab,'symbol'])
colnames(nodes) = paste('node',colnames(nodes),sep='.')
nodes=cbind(nodes, 
            angle=ifelse(posB[nodes$node.lab,'strand']=="+",90,-90), 
            radius=posB[nodes$node.lab,'end']-posB[nodes$node.lab,'start'], 
            offset=0, #posB[nodes$node.lab,'start']-posB[nodes$node.lab,'end'],
            hjust=ifelse(posB[nodes$node.lab,'strand']=="+",1,0),
            vjust=ifelse(posB[nodes$node.lab,'strand']=="+",0,1),
            # extra column
            strand=posB[nodes$node.lab,'strand'],
            type=posB[nodes$node.lab,'type']
)
# only show protein-coding
#nodes = nodes[posB[nodes$node.lab,'type']=='protein_coding',]
#nodes = nodes[posB[nodes$node.lab,'symbol'] %in% c('MAPT','KANSL1'),]
write.csv(nodes,"/tmp/nodes.csv",quote=F, row.names = F)

pdf(paste(gsub("[:-]","_",gsub(",","",regions)), "hiveplot","pdf", sep="."), width = 10, height=10)
# background in white
hive1$nodes$color[grep("ENSG", hive1$nodes$lab)]='black'
plotHive(hive1, ch=(ranger[2]-ranger[1])/10, method='abs', 
         axLabs = c("enhancers", "genes"),
         axLab.pos = as.integer(c((ranger[2]-ranger[1])/8,(ranger[2]-ranger[1])/8)),
         axLab.gpar = gpar(col = c("green","blue"), fontsize = 10),
         bkgnd = "NA",
         anNodes = "/tmp/nodes.csv",
         anNode.gpar = gpar(col = ifelse(nodes$type=='protein_coding','#000000','#00000066'), fontsize=5)
)
# background in black
hive1$nodes$color[grep("ENSG", hive1$nodes$lab)]='white'
plotHive(hive1, ch=(ranger[2]-ranger[1])/10, method='abs', 
         axLabs = c("enhancers", "genes"),
         axLab.pos = as.integer(c((ranger[2]-ranger[1])/8,(ranger[2]-ranger[1])/8)),
         axLab.gpar = gpar(col = c("green","blue"), fontsize = 10),
         bkgnd = "black",
         anNodes = "/tmp/nodes.csv",
         anNode.gpar = gpar(col = ifelse(nodes$type=='protein_coding','#ffffff','#ffffff66'), fontsize=5)
)
dev.off()