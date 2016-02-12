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
if (!require("HiveR")) {
  install.packages("HiveR", dependencies = TRUE)
  library(HiveR)
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
regions="chr17:43,583,680-44,506,585"; bins=500;

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
colnames(posB) = c('chr','start','end','name','score','strand')
rownames(posA) = posA$name; rownames(posB) = posB$name

region = strsplit(gsub(",","",regions), "[:-]")[[1]]
c=region[1]; s=as.numeric(region[2]); e=as.numeric(region[3])
step=(e-s)/bins

message("calculating the gene-enhancer correlation matrix ...")
# ===============================================
nmsA = subset(posA, chr==c & start<e & end>s, select=name)[,1]
nmsB = subset(posB, chr==c & start<e & end>s, select=name)[,1]
df = cor(t(valA[nmsA,]), t(valB[nmsB,]), method='spearman') # row is A and col is B

CUTOFF=0.4

# filter abs(cor())<0.5
df[is.na(df)]=0
df[abs(df)<CUTOFF]=0

message("plot the correlation matrix ...")
# ===============================================
# gene-gene interaction
hive1 <- adj2HPD(df, axis.cols =c('gray','gray'), desc = "Gene-enhancer correaltion")  # now hive1$nodes are first rows and then cols
# set axis
hive1$nodes$axis=as.integer(c(rep(1,nrow(df)),rep(2,ncol(df))))
# set radius of nodes
hive1$nodes$radius=c(ifelse(posA[nmsA, 'strand']=="-",posA[nmsA,'end'],posA[nmsA,'start']), ifelse(posB[nmsB, 'strand']=="-",posB[nmsB,'end'],posB[nmsB,'start']))
hive1$nodes$radius=hive1$nodes$radius-min(hive1$nodes$radius)+100
# set size of nodes
hive1$nodes$size=.5
# set color of nodes
hive1$nodes$color=c(rep('magenta',nrow(df)), rep('black',ncol(df)))
# set color of edge
hive1$edges$color=ifelse(hive1$edges$weight>0, rgb(1,0,0,abs(hive1$edges$weight)), rgb(1,1,1,abs(hive1$edges$weight)))

## add two artificial nodes on each axis, min-100 and max+100
ranger=range(hive1$nodes$radius);
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
hive1$edges$weight=.5

# remove gene self correlation (e.g. R=1)
hive1 = mineHPD(hive1, option = "remove zero edge")

pdf(paste(gsub(",","",regions), "hiveplot","pdf", sep="."), width = 4, height=4)
chkHPD(hive1)
plotHive(hive1, ch=.1, method='norm', 
         axLabs = c("enhancers", "genes"),
         bkgnd = "white",
         anNodes = "/tmp/nodes.txt",
         anNode.gpar = gpar(col = "black")
)


node.lab,node.text,angle,radius,offset,hjust,vjust
eye_33, "green", 240, 2.5, 0.2, 0.5, 0.5
eye_37, "hazel", 240, 2.5, 0.2, 0.5, 0.5
eye_41, "blue", 240, 2.5, 0.2, 0.5, 0.5
eye_45, "brown", 240, 2.5, 0.2, 0.5, 0.5
hair_1, "black", 60, 2.0, 0.2, 0.5, 0.5
hair_2, "brown", 60, 2.0, 0.2, 0.5, 0.5
hair_3, "red", 60, 2.0, 0.2, 0.5, 0.5
hair_4, "blond", 60, 2.0, 0.2, 0.5, 0.5


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
image(1:bins, 1:bins, df, zlim=c(-1,1))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

# creates a 5 x 5 inch image
png(paste(gsub(",","",regions),"correlationplot","png", sep=".")),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

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


