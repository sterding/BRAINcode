## HTNE vs. gene: which one is better for cell-type specificity? e.g. CELF2, MAPT locus



valA=args[1]  # e.g. eRNA.expression.tab
valB=args[2]  # e.g. genes.expression.tab
posA=args[3]  # e.g. eRNA postion bed
posB=args[4]  # e.g. genes postition bed
regions=args[5]
CUTOFF=args[6]

if(is.na(CUTOFF)) CUTOFF=0.1

# debug
valA="eRNA.meanRPM.xls"; valB="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls"; 
posA="eRNA.bed"; posB="~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed";
regions="chr17_43583680_44506585";
covURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"

message("reading data ...")
# ===============================================
valA=read.table(valA, header=T, stringsAsFactors =F, check.names =F)  # e.g. eRNA.expression.tab
rownames(valA) = valA[,1]; valA=valA[,-1]; 

valB=read.table(valB, header=T, stringsAsFactors =F, check.names =F)  # e.g. genes.expression.tab
rownames(valB) = valB[,1]; valB=valB[,-1]; 
colnames(valB)=gsub("FPKM.","",colnames(valB))
colnames(valB)=gsub("_0$","",colnames(valB))

posB=read.table(posB, header=F, stringsAsFactors =F)  # e.g. genes postition bed
colnames(posB) = c('chr','start','end','name','score','strand','symbol','type')
rownames(posB) = posB$name

posA=read.table(posA, header=F, stringsAsFactors =F)  # e.g. eRNA postion bed
colnames(posA) = c('chr','start','end','name') #,'score','strand')
posA$score=0; posA$strand=".";
rownames(posA) = posA$name;

message("cleanup data ...")
# ===============================================
# remove controls/outliers etc.
require(RCurl)
cov=read.delim(textConnection(getURL(covURL)), stringsAsFactors =F)
head(cov)
selected = cov$sampleName[cov$BRAINCODE.final.selection==1]
colnames(valB)
valB = valB[,selected]

# convert to long format
valB.long = valB
require('reshape2')
# convert sample to catelog name, e.g. HC_ND34770_FB_6_rep1 --> HC_FB
colnames(valB.long) = gsub("(.*)_.*(_.*)_[0-9]_rep.*","\\1\\2",colnames(valB.long))
# HC_SNDA and ILB_SNDA --> HCILB_SNDA
colnames(valB.long)[colnames(valB.long) %in% c('ILB_SNDA','HC_SNDA')]="HCILB_SNDA"
valB.long = melt(as.matrix(valB.long)) # melting a matix to include rownames: http://stackoverflow.com/questions/19826832/why-reshape2s-melt-cannot-capture-rownames-in-the-transformation
head(valB.long)
colnames(valB.long) = c("id","group",'fpkm')
## convert to log scale to stablize variance
valB.long$log10fpkm = log10(valB.long$fpkm + 0.01)

# calcluate mean, se, etc. per group per gene
library(plyr)
valB.stat = ddply(valB.long, c('id','group'), summarise,
      N = length(fpkm),
      mean = mean(fpkm),
      sd = sd(fpkm),
      se = sd / sqrt(N)
)
# rename the levels of group
valB.stat$group = factor(valB.stat$group, levels=c('HCILB_SNDA','HC_TCPY','HC_MCPY','HC_FB','HC_PBMC'))
head(valB.stat$group)
# add gene symbol 
valB.stat$symbol = posB$symbol[match(valB.stat$id, posB$name)]
head(valB.stat)

# plot some genes
df=subset(valB.stat, symbol %in% c('LINC00710','CELF2', 'MAPT','KANSL1', 'DCBLD2','NR4A2','TH'))
df=subset(valB.stat, symbol %in% c('CELF2','LINC00710','CELF2-AS2'))
df=subset(valB.stat, symbol %in% c('MAPT','CELF2','KANSL1','LINC00710'))

df
library(ggplot2)
p <- ggplot(df, aes(colour=symbol, y=mean, x=group))
p + geom_line(aes(group=symbol)) + geom_point()+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + ylab("FPKM") + theme_bw()
#p + geom_line(aes(group=symbol)) + geom_point()+ geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, colour='black') + scale_y_log10() + ylab("log10(FPKM)") + theme_bw()

# same order of columns for A and B
common = intersect(colnames(valA), colnames(valB))
valA=valA[,common]; valB=valB[,common]



region = strsplit(gsub(",","",regions), "[:-_]")[[1]]
c=region[1]; s=as.numeric(region[2]); e=as.numeric(region[3])
