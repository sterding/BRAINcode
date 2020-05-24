# script to plot correlation/interaction plot

args<-commandArgs(TRUE)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

valA=args[1]  # e.g. eRNA.expression.tab
valB=args[2]  # e.g. genes.expression.tab
posA=args[3]  # e.g. eRNA postion bed
posB=args[4]  # e.g. genes postition bed
regions=args[5]
bins=args[6] # 100

# debug
valA="eRNA.meanRPM.xls"; valB="/data/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls"; 
posA="eRNA.bed"; posB="/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed";
regions="chr17_43583680_44506585"; bins=500;

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

region = strsplit(gsub(",","",regions), "[:-_]")[[1]]
c=region[1]; s=as.numeric(region[2]); e=as.numeric(region[3])
step=(e-s)/bins

message("calculating the binned correlation matrix ...")
# ===============================================
df=matrix(data=0, nrow=bins, ncol=bins)
for(i in 1:bins){
  Si=start+(i-1)*step; Ei=Si+step;
  nms = subset(posA, chr==chr & start<Ei & end>Si, select=name)[,1]
  X=apply(valA[nms,],2,mean)
  message(paste("   ",i,"bin..."))
  for(j in 1:bins){
    Sj=start+(j-1)*step; Ej=Sj+step;
    nms = subset(posB, chr==chr & start<Ej & end>Sj, select=name)[,1]
    Y=apply(valB[nms,],2,mean)
    df[i,j]=cor(X,Y,method='spearman')
  }
}

message("plot the correlation matrix ...")
# ===============================================

# creates a 5 x 5 inch image
png(paste(gsub(",","",regions),"correlationplot","png", sep="."),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

heatmap.2(df,
          main = "Correlation", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          Colv = F,
          Rowv = F,
          symbreaks=T)
dev.off()