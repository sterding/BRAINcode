## HTNE vs. gene: which one is better for cell-type specificity? e.g. CELF2, MAPT locus

# debug
regions="chr17_43583680_44506585";
gene='CELF2'

expGene="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls"; 
posGene="~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed";
covURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"

message("reading data ...")
# ===============================================
expGene=read.table(expGene, header=T, stringsAsFactors =F, check.names =F)  # e.g. genes.expression.tab
rownames(expGene) = expGene[,1]; expGene=expGene[,-1]; 
colnames(expGene)=gsub("FPKM.","",colnames(expGene))
colnames(expGene)=gsub("_0$","",colnames(expGene))

posGene=read.table(posGene, header=F, stringsAsFactors =F)  # e.g. genes postition bed
colnames(posGene) = c('chr','start','end','name','score','strand','symbol','type')
rownames(posGene) = posGene$name

message("cleanup data ...")
# ===============================================
# colors
require(RCurl)
cl=read.delim(textConnection(getURL("https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/pub?gid=1995457670&single=true&output=tsv")), stringsAsFactors =F)
head(cl)
colours=paste0("#",cl$HEX[match(c('HCILB_SNDA','HC_PY', 'HC_nonNeuron', 'HC_TCPY','HC_MCPY','HC_FB','HC_PBMC'), cl$ITEM)])
names(colours) = c('HCILB_SNDA','HC_PY', 'HC_nonNeuron', 'HC_TCPY','HC_MCPY','HC_FB','HC_PBMC')

# remove controls/outliers etc.
require(RCurl)
cov=read.delim(textConnection(getURL(covURL)), stringsAsFactors =F)
head(cov)
selected = cov$sampleName[cov$BRAINCODE.final.selection==1]
colnames(expGene)
expGene = expGene[,selected]

# convert to long format
expGene.long = expGene
require('reshape2')
# convert sample to catelog name, e.g. HC_ND34770_FB_6_rep1 --> HC_FB
colnames(expGene.long) = gsub("(.*)_.*(_.*)_[0-9]_rep.*","\\1\\2",colnames(expGene.long))
# HC_SNDA and ILB_SNDA --> HCILB_SNDA
colnames(expGene.long)[colnames(expGene.long) %in% c('ILB_SNDA','HC_SNDA')]="HCILB_SNDA"
expGene.long = melt(as.matrix(expGene.long)) # melting a matix to include rownames: http://stackoverflow.com/questions/19826832/why-reshape2s-melt-cannot-capture-rownames-in-the-transformation
head(expGene.long)
colnames(expGene.long) = c("id","group",'fpkm')
## convert to log scale to stablize variance
expGene.long$log10fpkm = log10(expGene.long$fpkm + 0.01)

# major group
head(expGene.long)
expGene.long$group3=ifelse(grepl("SNDA",expGene.long$group),'HCILB_SNDA',ifelse(grepl("PY",expGene.long$group),"HC_PY","HC_nonNeuron"))

# calcluate mean, se, etc. per group per gene
library(plyr)
expGene.stat = ddply(expGene.long, c('id','group'), summarise,
      N = length(fpkm),
      mean = mean(fpkm),
      sd = sd(fpkm),
      se = sd / sqrt(N)
)
# rename the levels of group
expGene.stat$group = factor(expGene.stat$group, levels=c('HCILB_SNDA','HC_TCPY','HC_MCPY','HC_FB','HC_PBMC'))
# add gene symbol 
expGene.stat$symbol = posGene$symbol[match(expGene.stat$id, posGene$name)]

expGene.stat3 = ddply(expGene.long, c('id','group3'), summarise,
                     N = length(fpkm),
                     mean = mean(fpkm),
                     sd = sd(fpkm),
                     se = sd / sqrt(N)
)
# rename the levels of group
expGene.stat3$group3 = factor(expGene.stat3$group3, levels=c('HCILB_SNDA','HC_PY','HC_nonNeuron'))
# add gene symbol 
expGene.stat3$symbol = posGene$symbol[match(expGene.stat3$id, posGene$name)]

# plot some genes
df=subset(expGene.stat, symbol %in% c('LINC00710','CELF2', 'MAPT','KANSL1', 'DCBLD2','NR4A2','TH'))
df=subset(expGene.stat, symbol %in% c('CELF2','LINC00710','CELF2-AS2'))
df=subset(expGene.stat, symbol %in% c('MAPT','CELF2','KANSL1','LINC00710'))

pdf(paste0("~/eRNAseq/eRNA.vs.gene.",gene,".pdf"), width=3, height=4)
library(ggplot2)
df=subset(expGene.stat, symbol %in% gene)
p <- ggplot(df, aes(colour=symbol, y=mean, x=group, fill=group))
#p + geom_line(aes(group=symbol)) + geom_point(shape=16,size=3, fill="red")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("gene FPKM")+theme_bw() + theme(legend.position="top") 
# tips for assigning color manually: http://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2
p + geom_bar(stat="identity", color=NA) + scale_fill_manual(values=colours) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("gene FPKM") + theme_bw() + theme_classic()  + theme(legend.position="none") 

df=subset(expGene.stat3, symbol %in% gene)
p <- ggplot(df, aes(colour=symbol, y=mean, x=group3, fill=group3))
#p + geom_line(aes(group=symbol)) + geom_point(shape=16,size=3, fill="red")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("gene FPKM")+theme_bw() +theme(legend.position="top") 
p + geom_bar(stat="identity", color=NA) + scale_fill_manual(values=colours) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("gene FPKM") + theme_bw() + theme_classic() + theme(legend.position="none") 

## intron meanRPM of genes
## ======================================
metaIntron="~/neurogen/rnaseq_PD/results/merged/metaIntron.meanRPM.allSamples.xls"
metaIntron=read.table(metaIntron, header=T, stringsAsFactors =F, check.names =F)  # e.g. genes.expression.tab
rownames(metaIntron) = metaIntron[,1]; metaIntron=metaIntron[,-1]; 
metaIntron = metaIntron[,selected]
metaIntron.long = metaIntron
require('reshape2')
colnames(metaIntron.long) = gsub("(.*)_.*(_.*)_[0-9]_rep.*","\\1\\2",colnames(metaIntron.long))
colnames(metaIntron.long)[colnames(metaIntron.long) %in% c('ILB_SNDA','HC_SNDA')]="HCILB_SNDA"
metaIntron.long = melt(as.matrix(metaIntron.long)) # melting a matix to include rownames: http://stackoverflow.com/questions/19826832/why-reshape2s-melt-cannot-capture-rownames-in-the-transformation
head(metaIntron.long)
colnames(metaIntron.long) = c("id","group",'fpkm')
metaIntron.long$log10fpkm = log10(metaIntron.long$fpkm + 0.001)
metaIntron.long$group3=ifelse(grepl("SNDA",metaIntron.long$group),'SNDA',ifelse(grepl("PY",metaIntron.long$group),"PY","Non-neuron"))

library(plyr)
metaIntron.stat = ddply(metaIntron.long, c('id','group'), summarise,
                     N = length(fpkm),
                     mean = mean(fpkm),
                     sd = sd(fpkm),
                     se = sd / sqrt(N)
)
# rename the levels of group
metaIntron.stat$group = factor(metaIntron.stat$group, levels=c('HCILB_SNDA','HC_TCPY','HC_MCPY','HC_FB','HC_PBMC'))
# add gene symbol 
metaIntron.stat$symbol = posGene$symbol[match(metaIntron.stat$id, posGene$name)]

metaIntron.stat3 = ddply(metaIntron.long, c('id','group3'), summarise,
                        N = length(fpkm),
                        mean = mean(fpkm),
                        sd = sd(fpkm),
                        se = sd / sqrt(N)
)
# rename the levels of group
metaIntron.stat3$group3 = factor(metaIntron.stat3$group3, levels=c('SNDA','PY','Non-neuron'))
# add gene symbol 
metaIntron.stat3$symbol = posGene$symbol[match(metaIntron.stat3$id, posGene$name)]

# plot some genes
df=subset(metaIntron.stat, symbol %in% gene)


library(ggplot2)
p <- ggplot(df, aes(colour=symbol, y=mean, x=group))
p + geom_line(aes(group=symbol)) + geom_point(shape=16,size=3, fill="red")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("Intronic meanPKM") + theme_bw() + theme(legend.position="top")

df=subset(metaIntron.stat3, symbol %in% gene)
p <- ggplot(df, aes(colour=symbol, y=mean, x=group3))
p + geom_line(aes(group=symbol)) + geom_point(shape=16,size=3, fill="red")+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("Intronic meanPKM") + theme_bw() + theme(legend.position="top") 

## meanRPM per group
metaIntron=c()
celltypes=c("HCILB_SNDA",'HC_TCPY','HC_MCPY','HC_FB','HC_PBMC')
for(i in celltypes){
  df=read.table(paste0("~/eRNAseq/",i,"/meanRPM.of.metaintron.by.gene.tab"), header=T, stringsAsFactors =F);
  df$symbol = posGene$symbol[match(df$id, posGene$name)]
  metaIntron = c(metaIntron, df$meanRPM[df$symbol==gene])
}
names(metaIntron) = celltypes
barplot(metaIntron, col=colours[celltypes],cex.names=0.5, ylab="Normalized expression of intron")

# major
metaIntron=c();
celltypes=c("HCILB_SNDA",'HC_PY','HC_nonNeuron')
for(i in celltypes){
  df=read.table(paste0("~/eRNAseq/",i,"/meanRPM.of.metaintron.by.gene.tab"), header=T, stringsAsFactors =F);
  df$symbol = posGene$symbol[match(df$id, posGene$name)]
  metaIntron = c(metaIntron, df$meanRPM[df$symbol==gene])
}
names(metaIntron) = celltypes
barplot(metaIntron, col=colours[celltypes],cex.names=0.5, ylab="Normalized expression of intron")

## length of the longest HTNE cluster
## ======================================
clusterHTNE=c()
celltypes = c("HCILB_SNDA",'HC_TCPY','HC_MCPY','HC_FB','HC_PBMC') # minor
for(i in celltypes){
  df=read.table(paste0("~/eRNAseq/",i,"/eRNA.cluster.bed"), header=F, stringsAsFactors =F);
  clusterHTNE = c(clusterHTNE, max(df$V5[grep(gene, df$V6,fixed=T)]))
}
clusterHTNE=clusterHTNE/1000
names(clusterHTNE) = celltypes
#plot(clusterHTNE, type='b',xaxt="n",ylab="Length of HTNE cluster (Kb)")
#points(1:length(clusterHTNE),clusterHTNE,pch=19,col=colours[celltypes])
#axis(1, at=1:length(clusterHTNE), cex.axis=0.7, labels=celltypes)
barplot(clusterHTNE, col=colours[celltypes],cex.names=0.5, ylab="Length of HTNE cluster (Kb)")

# majpr
clusterHTNE=c()
celltypes = c("HCILB_SNDA",'HC_PY','HC_nonNeuron') # major
for(i in celltypes){
  df=read.table(paste0("~/eRNAseq/",i,"/eRNA.cluster.bed"), header=F, stringsAsFactors =F);
  clusterHTNE = c(clusterHTNE, max(df$V5[grep(gene, df$V6,fixed=T)]))
}
clusterHTNE=clusterHTNE/1000
names(clusterHTNE) = celltypes
#plot(clusterHTNE, type='b',xaxt="n", ylab="Length of HTNE cluster (Kb)")
#points(1:length(clusterHTNE),clusterHTNE,pch=19,col=colours[celltypes])
#axis(1, at=1:length(clusterHTNE), cex.axis=0.7, labels=celltypes)
barplot(clusterHTNE, col=colours[celltypes],cex.names=0.5, ylab="Length of HTNE cluster (Kb)")


dev.off()
