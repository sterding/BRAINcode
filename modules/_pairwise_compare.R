###########################################
# Rscript to check the consistency of RNAseq between two samples
# Usage: Rscript $pipeline_path/modules/_pairwise_compare.R HC_M0235-4_PBMC_6_rep1.amplified HC_M0235-4_PBMC_6_rep1.unamplified
# Author: Xianjun Dong
# Version: 0.0
# Date: 2015-Apr-12
###########################################
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(ggplot2, ggExtra)
library(ggplot2)
library(ggExtra) # install.packages('ggExtra')

args<-commandArgs(TRUE)
sample1=args[1]  
sample2=args[2]
FPKMfile = ifelse(is.na(args[3]), "~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls", args[3])

addmarkers=0
filter=0; 
PSEUDOCOUNT=1e-2

# test:
# sample1='HC_UWA616_SN_6_rep1.amplified'; sample2='HC_UWA616_SN_6_rep1.unamplified'; FPKMfile ="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls"; addmarkers=0; PSEUDOCOUNT=1e-5;filter=0.001

if(grepl("_SN_",sample1) & grepl("_SNDA_",sample2)) addmarkers=1

# setwd("~/neurogen/rnaseq_PD/results/merged/"); sample1="HC_UWA616_SN_6_rep1.amplified"; sample2="HC_UWA616_SNDA_2_rep1"; FPKMfile="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls";PSEUDOCOUNT=1e-2
print(paste(sample1, sample2, FPKMfile, "PSEUDO:",PSEUDOCOUNT,"FILTER:",filter, "ADDMARKER:",addmarkers))
message("loading data...")

if(tools::file_ext(FPKMfile)=='rds') fpkm=readRDS(FPKMfile) else {
  fpkm=read.table(FPKMfile, header=T, row.names = 1, check.names =F, stringsAsFactors =F);  # table with header (1st row) and ID (1st column)
}
#fpkm=fpkm[,grep("FPKM",colnames(fpkm))]; 
colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))

if(sample1 %in% colnames(fpkm)) x=fpkm[,sample1] else {message(paste("cannot find the sample", sample1, "in the input expression file")); quit('no');}
if(sample2 %in% colnames(fpkm)) y=fpkm[,sample2] else {message(paste("cannot find the sample", sample2, "in the input expression file")); quit('no');}

isna=is.na(x) | is.na(y); x=x[!isna]; y=y[!isna]; 

if(filter>0){
  to_keep = x>=filter & y>=filter
  x=x[to_keep]; y=y[to_keep];
}

# marker genes
if(addmarkers==1){
    neuron=fpkm$gene_short_name %in% c('RBFOX3','DRD2','ENO3', 'TH','SYT1','SLC12A5', 'ENO2')
    microglia=fpkm$gene_short_name %in% c('P2RY12', 'GPR34','CSF1R','CD53')
    glia=fpkm$gene_short_name %in% c('ALDH1L1','MAG','GFAP','PMP22','OLIG2','OLIG1','AIF1','MOG','SOX10','GJC2','GJB6')
    selected_markers=neuron | microglia | glia
    cols=ifelse(selected_markers & neuron, 'red', ifelse(selected_markers & microglia, 'blue', ifelse(selected_markers & glia, 'lightblue', NA)))
    # length(cols)
    selected_markers=selected_markers[!isna]
    cols=cols[!isna & !is.na(cols)]
}

message("ploting data...")

x=log10(x + PSEUDOCOUNT); y=log10(y + PSEUDOCOUNT);

pdf(paste("xyplot", sample1, "vs", sample2, "PSEUDO",PSEUDOCOUNT,"FILTER",filter, "pdf", sep="."), width=5, height=5, paper='us', useDingbats=FALSE)
par(pty="s"); # to make sure the frame is a square (width=height)
plot(unique(round(cbind(x,y),3)), pch='.', cex=1, main=basename(FPKMfile), xlab=paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")"), ylab=paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"), xlim=range(x,y), ylim=range(x,y), xaxs="r", yaxs="r")
if(addmarkers==1){
    points(x=x[selected_markers], y=y[selected_markers], col=cols, pch=19, cex=1)
    text(x=x[selected_markers],y=y[selected_markers],labels = fpkm$gene_short_name[selected_markers], adj=c(0,1), col=cols)
}
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(x,y), 3)), bty='n')

# heatmap
smoothScatter(x,y,xlab=paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")"), ylab=paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"), xlim=range(x,y), ylim=range(x,y), xaxs="r", yaxs="r")
if(addmarkers==1){
  points(x=x[selected_markers], y=y[selected_markers], col=cols, pch=19, cex=1)
  text(x=x[selected_markers],y=y[selected_markers],labels = fpkm$gene_short_name[selected_markers], adj=c(0,1), col=cols)
}
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(x,y), 3)), bty='n')

# add density curve
df=data.frame(x=x,y=y)
p=ggplot(df, aes(x,y))+ geom_point(size=.6, alpha=0.6) + xlab(paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")")) + ylab(paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"))
p=ggMarginal(p)
print(p)

dev.off()

quit('no')



## SN vs. SNDA
setwd("~/neurogen/rnaseq_PD/results/merged/"); sample1="HC_UWA616_SN_6_rep1.amplified"; sample2="HC_UWA616_SNDA_2_rep1"; 
FPKMfile="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cufflinks.allSamples.uniq.xls";PSEUDOCOUNT=1e-2
fpkm=read.table(FPKMfile, header=T, check.names =F, stringsAsFactors =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))
if(sample1 %in% colnames(fpkm)) x=fpkm[,sample1] else {message(paste("cannot find the sample", sample1, "in the input expression file")); quit('no');}
if(sample2 %in% colnames(fpkm)) y=fpkm[,sample2] else {message(paste("cannot find the sample", sample2, "in the input expression file")); quit('no');}
isna=is.na(x) | is.na(y); x=x[!isna]; y=y[!isna]; 

message("ploting data...")

pdf(paste("xyplot", sample1, "vs", sample2, PSEUDOCOUNT, "pdf", sep="."), width=5, height=5, paper='us', useDingbats=FALSE)
par(pty="s"); # to make sure the frame is a square (width=height)
# marker genes
markers=list(
    neuron_general=c('RBFOX3','MAP2','TUBB3','ENO3','SYT1','SLC12A5', 'ENO2','NEFL'), # RBFOX3=NeuN; SYT1=Synaptotagmin 1; SLC12A5=KCC2; NEFL=Neurofilament
    neuron_dopamine=c('DRD2','TH','NR4A2','DDC','SLC6A3','SLC18A1'), # NR4A2=NURR1
    glia_micro=c('P2RY12', 'GPR34','CSF1R','CD53','HEXB','ITGAM','PTPRC'), # ITGAM=CD11b; PTPRC=CD45
    glia_astro=c('ALDH1L1','GFAP','GJB6'), # GJB6=Connexin 30;
    glia_oligo=c('MAG','OLIG2','OLIG1','SOX10','GJC2'), # GJC2=connexin-47
    glia_myeli=c('PMP22','MOG','MPZ','MBP')
)
selected_markers=rev(match(unlist(markers), fpkm$gene_short_name))
col=c('darkred','magenta','forestgreen','blue','cadetblue1','cornflowerblue')
cols=rev(rep(col, sapply(markers, length)))
X=log10(x + PSEUDOCOUNT); Y=log10(y + PSEUDOCOUNT);

plot(unique(round(cbind(X,Y),3)), pch='.', cex=0.6, main=basename(FPKMfile), xlab=paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")"), ylab=paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"), xlim=range(X,Y), ylim=range(X,Y), xaxs="r", yaxs="r")
points(x=X[selected_markers], y=Y[selected_markers], col=cols, pch=19, cex=1)
text(x=X[selected_markers],y=Y[selected_markers],labels = fpkm$gene_short_name[selected_markers], adj=c(0,1), col=cols)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(X,Y), 3)), bty='n')
legend("bottomright", names(markers), col=col, pch=19, cex=0.6, bty='n')

# barplot of fold-change
X=log2(x + PSEUDOCOUNT); Y=log2(y + PSEUDOCOUNT);
barplot(Y[selected_markers]-X[selected_markers], xlab="log2(fold-change of SNDA/SN)",horiz =T, col=cols, names.arg = fpkm$gene_short_name[selected_markers], cex.names =0.5, las=1, border=NA)
legend("topright", names(markers), col=col, pch=15, cex=0.6, bty='n')

# simplified version
markers=list(
    #neuron_general=c('RBFOX3','ENO3','SLC12A5', 'ENO2'), # RBFOX3=NeuN; SYT1=Synaptotagmin 1; SLC12A5=KCC2; NEFL=Neurofilament
    neuron_dopamine=c('TH','SLC6A3','SLC18A1','DRD2'), # NR4A2=NURR1; SLC6A3=DAT; SLC18A1=VMAT1
    glia_micro=c('P2RY12','PTPRC'), # ITGAM=CD11b; PTPRC=CD45
    glia_astro=c('GFAP','GJB6'), # GJB6=Connexin 30;
    glia_oligo=c('OLIG2','OLIG1'), # GJC2=connexin-47
    glia_myeli=c('PMP22','MOG')
)
selected_markers=rev(match(unlist(markers), fpkm$gene_short_name))
col=c('magenta','forestgreen','blue','cadetblue1','cornflowerblue')
cols=rev(rep(col, sapply(markers, length)))
X=log10(x + PSEUDOCOUNT); Y=log10(y + PSEUDOCOUNT);

message("saving data...")
write.table(fpkm[selected_markers,c('gene_short_name',sample1,sample2)],file=paste("xyplot", sample1, "vs", sample2, "RPKM", "xls", sep="."), sep="\t", col.names = T)

plot(unique(round(cbind(X,Y),3)), pch='.', cex=0.6, main=basename(FPKMfile), xlab=paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")"), ylab=paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"), xlim=range(X,Y), ylim=range(X,Y), xaxs="r", yaxs="r")
points(x=X[selected_markers], y=Y[selected_markers], col=cols, pch=19, cex=1)
text(x=X[selected_markers],y=Y[selected_markers],labels = fpkm$gene_short_name[selected_markers], adj=c(0,1), col=cols)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(X,Y), 3)), bty='n')
legend("bottomright", names(markers), col=col, pch=19, cex=0.6, bty='n')

X=log2(x + PSEUDOCOUNT); Y=log2(y + PSEUDOCOUNT);
barplot(Y[selected_markers]-X[selected_markers], xlab="log2(fold-change of SNDA/SN)",horiz =T, col=cols, names.arg = fpkm$gene_short_name[selected_markers], cex.names =0.5, las=1, border=NA, axes=F, xlim=c(-6,3))
axis(1,at=c(-6:3),labels = c(-6:3))
legend("topleft", names(markers), col=col, pch=19, cex=0.8, bty='n')

dev.off()




## hocky curve investigation (for Clemens)

FPKMfile="/data/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls";
fpkm=read.table(FPKMfile, header=T, check.names =F, stringsAsFactors =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; fpkm=fpkm[,grep("FPKM",colnames(fpkm))]; colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))

amp1="HC_M0235-4_PBMC_6_rep1.amplified"; non1="HC_M0235-4_PBMC_6_rep1.unamplified";
amp2="HC_UWA616_SN_6_rep1.amplified"; non2="HC_UWA616_SN_6_rep1.unamplified";

if(amp1 %in% colnames(fpkm)) AMP1=fpkm[,amp1] else {message(paste("cannot find the sample", amp1, "in the input expression file")); quit('no');}
if(non1 %in% colnames(fpkm)) NON1=fpkm[,non1] else {message(paste("cannot find the sample", non1, "in the input expression file")); quit('no');}
if(amp2 %in% colnames(fpkm)) AMP2=fpkm[,amp2] else {message(paste("cannot find the sample", amp2, "in the input expression file")); quit('no');}
if(non2 %in% colnames(fpkm)) NON2=fpkm[,non2] else {message(paste("cannot find the sample", non2, "in the input expression file")); quit('no');}

isna=is.na(AMP1) | is.na(NON1) | is.na(AMP2) | is.na(NON2);
AMP1=AMP1[!isna]; NON1=NON1[!isna]; AMP2=AMP2[!isna]; NON2=NON2[!isna]

message("ploting data...")

pdf("xyplot.amp.vs.nonamp.pdf", width=6, height=6, paper='us')

rang=range(AMP1,NON1,AMP2,NON2)+PSEUDOCOUNT;

plot(AMP1 + PSEUDOCOUNT, NON1 + PSEUDOCOUNT, log='xy',pch='.', asp=1, cex=0.6, main="AMP1 vs. NON1", xlab=amp1, ylab=non1, xlim=rang, ylim=rang)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(log10(AMP1 + PSEUDOCOUNT), log10(NON1 + PSEUDOCOUNT)), 3)), bty='n')

plot(AMP2 + PSEUDOCOUNT, NON2 + PSEUDOCOUNT, log='xy',pch='.', asp=1, cex=0.6, main="AMP2 vs. NON2", xlab=amp2, ylab=non2, xlim=rang, ylim=rang)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(log10(AMP2 + PSEUDOCOUNT), log10(NON2 + PSEUDOCOUNT)), 3)), bty='n')

amp2="HC_B0254-4_PBMC_6_rep1"
if(amp2 %in% colnames(fpkm)) AMP2=fpkm[,amp2] else {message(paste("cannot find the sample", amp2, "in the input expression file")); quit('no');}
AMP2=AMP2[!is.na(AMP2)]
plot(AMP1 + PSEUDOCOUNT, AMP2 + PSEUDOCOUNT, log='xy',pch='.', asp=1, cex=0.6, main="AMP1 vs. AMP2", xlab=amp1, ylab=amp2, xlim=rang, ylim=rang)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(log10(AMP1 + PSEUDOCOUNT), log10(AMP2 + PSEUDOCOUNT)), 3)), bty='n')

amp1="HC_H1529-3_PBMC_6_rep1"; amp2="HC_H1560-2_PBMC_6_rep1"
if(amp1 %in% colnames(fpkm)) AMP1=fpkm[,amp1] else {message(paste("cannot find the sample", amp1, "in the input expression file")); quit('no');}
if(amp2 %in% colnames(fpkm)) AMP2=fpkm[,amp2] else {message(paste("cannot find the sample", amp2, "in the input expression file")); quit('no');}
AMP1=AMP1[!is.na(AMP1)]; AMP2=AMP2[!is.na(AMP2)]

plot(AMP1 + PSEUDOCOUNT, AMP2 + PSEUDOCOUNT, log='xy',pch='.', asp=1, cex=0.6, main="AMP1 vs. AMP2", xlab=amp1, ylab=amp2, xlim=rang, ylim=rang)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(log10(AMP1 + PSEUDOCOUNT), log10(AMP2 + PSEUDOCOUNT)), 3)), bty='n')

plot(NON1 + PSEUDOCOUNT, NON2 + PSEUDOCOUNT, log='xy',pch='.', asp=1, cex=0.6, main="NON1 vs. NON2", xlab=non1, ylab=non2, xlim=rang, ylim=rang)
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(log10(NON1 + PSEUDOCOUNT), log10(NON2 + PSEUDOCOUNT)), 3)), bty='n')

dev.off()
