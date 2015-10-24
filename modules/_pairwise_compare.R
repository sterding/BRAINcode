###########################################
# Rscript to check the consistency of RNAseq between two samples
# Usage: Rscript $pipeline_path/modules/_pairwise_compare.R HC_M0235-4_PBMC_6_rep1.amplified HC_M0235-4_PBMC_6_rep1.unamplified
# Author: Xianjun Dong
# Version: 0.0
# Date: 2015-Apr-12
###########################################

args<-commandArgs(TRUE)
sample1=args[1]  
sample2=args[2]
FPKMfile = ifelse(is.na(args[3]), "/data/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls", args[3])

PSEUDOCOUNT=1 #1e-5 

# setwd("~/neurogen/rnaseq_PD/results/merged/"); sample1="HC_M0235-4_PBMC_6_rep1.amplified"; sample2="HC_M0235-4_PBMC_6_rep1.unamplified"; FPKMfile="~/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls";PSEUDOCOUNT=1e-5
print(paste(sample1, sample2, FPKMfile))
message("loading data...")

fpkm=read.table(FPKMfile, header=T, check.names =F, stringsAsFactors =F);  # table with header (1st row) and ID (1st column)
rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]; fpkm=fpkm[,grep("FPKM",colnames(fpkm))]; colnames(fpkm)=gsub("FPKM.","",colnames(fpkm))

if(sample1 %in% colnames(fpkm)) x=fpkm[,sample1] else {message(paste("cannot find the sample", sample1, "in the input expression file")); quit('no');}
if(sample2 %in% colnames(fpkm)) y=fpkm[,sample2] else {message(paste("cannot find the sample", sample2, "in the input expression file")); quit('no');}

isna=is.na(x) | is.na(y); x=x[!isna]; y=y[!isna]

message("ploting data...")

x=log10(x + PSEUDOCOUNT); y=log10(y + PSEUDOCOUNT);

pdf(paste("xyplot", sample1, "vs", sample2, PSEUDOCOUNT, "pdf", sep="."), width=5, height=5, paper='us')
par(pty="s"); # to make sure the frame is a square (width=height)
plot(x, y, pch='.', cex=0.6, main=basename(FPKMfile), xlab=paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")"), ylab=paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"), xlim=range(x,y), ylim=range(x,y), xaxs="r", yaxs="r")
abline(a=0, b=1, col='red', lty=2, lwd=1)
legend("topleft", paste("Pearson's r =", round(cor(x,y), 3)), bty='n')
dev.off()

library("ggExtra"); library("ggplot2");
df=data.frame(x=x,y=y)
p=ggplot(df, aes(x,y))+ geom_point(size=.6, alpha=0.6) + xlab(paste(sample1, "log10(FPKM +",PSEUDOCOUNT,")")) + ylab(paste(sample2, "log10(FPKM +",PSEUDOCOUNT,")"))
p=ggMarginal(p)
ggsave(paste("xyplot+density", sample1, "vs", sample2, PSEUDOCOUNT, "pdf", sep="."), p)

quit('no')





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
