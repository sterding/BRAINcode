# Rscript for all coverage plots
# Update: Exclude the gap regions when calcluating genome size (2016.09.15)

setwd("~/neurogen/rnaseq_PD/results/coverage")

## ===============================================
## % of genome transcribed with different cutoff
## ===============================================

GENOMESIZE=3137161264  # hg19 genome size, see https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=&hgsid=512422197_TMKohgsX9wzqqDYa6ny5BLVqsYaD
GAP=239845127 # awk '{print $3-$2}' $GENOME/Annotation/Genes/hg19.gap.bed | datamash sum 1
genomesize=GENOMESIZE-GAP

df=read.table("coverageWithRPM.txt", header=T)

# optional: only show the 106 qualified samples
selected = readLines(pipe("cat ../merged/samplelist.HCILB_SNDA ../merged/samplelist.HC_PY ../merged/samplelist.HC_nonNeuron"))
df=subset(df, sample %in% selected)

df$celltype=gsub(".*_.*_(.*)_.*_.*","\\1",df$sample)
df$celltype[grepl("PY",df$celltype)]="PY"; df$celltype[grepl("PBMC|FB",df$celltype)]="NN"
df$sum=rowSums(df[,grep("RPM",colnames(df))])
df = df[with(df, order(celltype, -sum)), ]
df$sample <- factor(df$sample, unique(as.character(df$sample)))

# mean and sd per samples for RPM>=0.05
library(dplyr)
df %>% group_by(celltype) %>% mutate(RPM0.05 = RPMgt0.05 + RPMgt0.1 + RPMgt0.5 + RPMgt1) %>% summarise(mean=mean(RPM0.05*100/genomesize), median=median(RPM0.05*100/genomesize), sd=sd(RPM0.05*100/genomesize))

library(reshape2)
dflong=melt(df[,1:7], variable.name = "cutoff",value.name ="coverage")
#levels(dflong$cutoff)=rev(levels(dflong$cutoff))
library(ggplot2)
ggplot(dflong, aes(x=sample, y=100*coverage/genomesize, fill=cutoff, order = -as.numeric(cutoff))) +
    geom_bar(width=.5,position = position_stack(), stat="identity") +
    theme_bw() +
    ylab("Coverage of the whole genome (%)") +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=5), legend.justification=c(1,1), legend.position=c(1,1))
ggsave("coverageWithRPM.all.pdf", width=8, height=4)

df=read.table("coverageWithRPM.txt", header=T)
# optional: only show the 106 qualified samples
selected = readLines(pipe("cat ../merged/samplelist.HCILB_SNDA ../merged/samplelist.HC_PY ../merged/samplelist.HC_nonNeuron"))
df=subset(df, sample %in% selected)

df$celltype=gsub(".*_.*_(.*)_.*_.*","\\1",df$sample)
df$celltype[grepl("PY",df$celltype)]="PY"; df$celltype[grepl("PBMC|FB",df$celltype)]="NN"
df$sum=rowSums(df[,c("RPMgt1","RPMgt0.5","RPMgt0.1","RPMgt0.05")])
df = df[with(df, order(celltype, -sum)), ]
df$sample <- factor(df$sample, unique(as.character(df$sample)))

dflong=melt(df[,c("sample","RPMgt1","RPMgt0.5","RPMgt0.1","RPMgt0.05")], variable.name = "cutoff",value.name ="coverage")
#levels(dflong$cutoff)=rev(levels(dflong$cutoff))
ggplot(dflong, aes(x=sample, y=100*coverage/genomesize, fill=cutoff, order = -as.numeric(cutoff))) +
    geom_bar(width=.5,position = position_stack(), stat="identity") +
    theme_bw() +
    ylab("Coverage of the whole genome (%)") +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=5), legend.justification=c(1,1), legend.position=c(1,1))
ggsave("coverageWithRPM.RPMpt05.pdf", width=8, height=4)

## ===============================================
## three plots:
## 1. cummulative coverage of genome transcribed with 0.5 RPM
## 2. coverage of genome transcribed with 0.5 RPM for random 7 samples
## 3. cummulative coverage of genome transcribed with 5 reads
## ===============================================

pdf("coverage.cummulative.pdf", paper='us',width=5, height=4)

# 3 cell types in one plot
df=read.table("covered.0.05RPM.HCILB_SNDA.txt", header=F)
colnames(df)=c("samplecount", "coveredbp")
df=cbind(df, cumsum=cumsum(df$coveredbp), type='HCILB_SNDA')
par(mar=c(4,4,2,4), pty="s")
plot(df$cumsum, type='l', lwd=2, col='#F22A7B', ylim=c(0,2000000000), xlim=c(1,90), main="covered.0.05RPM", yaxt="n", xaxt='n', xlab='', ylab='')
df=read.table("covered.0.05RPM.HC_PY.txt", header=F)
colnames(df)=c("samplecount", "coveredbp")
df=cbind(df, cumsum=cumsum(df$coveredbp), type='HC_PY')
points(df$cumsum, type='l', lwd=2,col='#3182bd')
df=read.table("covered.0.05RPM.HC_nonNeuron.txt", header=F)
colnames(df)=c("samplecount", "coveredbp")
df=cbind(df, cumsum=cumsum(df$coveredbp), type='HC_nonNeuron')
points(df$cumsum, type='l', lwd=2,col='#513931')
axis(1, c(1,seq(20,90,20)), c(1,seq(20,90,20)), tck=0.02, mgp=c(3,0.2,0))
#axis(2, seq(0,2,.5)*1e9, labels=format(seq(0,2,.5),2), las=2, tck=0.02, mgp=c(3,0.2,0))
axis(2, genomesize*seq(0,60,10)/100, labels=round(genomesize*seq(0,60,10)/100/1e9,2), las=2)
axis(4, genomesize*seq(0,60,10)/100, labels=paste0(seq(0,60,10)), las=2)
abline(h=genomesize*seq(0,60,10)/100, lwd=.5, col='gray')
mtext("Covered percentage (%)", 4, line=2)
mtext("Sample count", 1, line=2)
mtext("Covered base pairs (in billion)", 2, line=3)

# coverage with same number of samples (N=5)
a=read.table("covered.0.05RPM.HCILB_SNDA.random5.txt", header=F)
a=100*apply(a,1,sum)/genomesize
df=data.frame(percentage=a,type='SNDA')
a=read.table("covered.0.05RPM.HC_PY.random5.txt", header=F)
a=100*apply(a,1,sum)/genomesize
df=rbind(df, data.frame(percentage=a,type='PY'))
a=read.table("covered.0.05RPM.HC_nonNeuron.random5.txt", header=F)
a=100*apply(a,1,sum)/genomesize
df=rbind(df, data.frame(percentage=a,type='NN'))
boxplot(percentage~type, df, outline =F, col=c('#F22A7B','#3182bd','#513931'), ylab="Covered percentage (%)", main="random 5 samples",yaxt="n")
axis(2, at=seq(10,30,5),labels=seq(10,30,5))
t.test(df$percentage[df$type=='SNDA'], df$percentage[df$type=='NN'])
t.test(df$percentage[df$type=='PY'], df$percentage[df$type=='NN'])
t.test(df$percentage[df$type=='SNDA'], df$percentage[df$type=='PY'])
library('dplyr')
df %>% group_by(type) %>% summarise(mean=mean(percentage), median=median(percentage), sd=sd(percentage))

# coverage with same number of samples (N=7)
a=read.table("covered.0.05RPM.HCILB_SNDA.random7.txt", header=F)
a=100*apply(a,1,sum)/genomesize
#a=try(system("cut -f2 random.covered.0.05RPM.HCILB_SNDA.*.txt | paste - - - - - - -",intern = T))
#a=100*apply(matrix(as.numeric(do.call(rbind, strsplit(a,"\t"))),nrow=100),1,sum)/genomesize
df=data.frame(percentage=a,type='SNDA')
#a=try(system("cut -f2 random.covered.0.05RPM.HC_PY.*.txt | paste - - - - - - -",intern = T))
#a=100*apply(matrix(as.numeric(do.call(rbind, strsplit(a,"\t"))),nrow=100),1,sum)/genomesize
a=read.table("covered.0.05RPM.HC_PY.random7.txt", header=F)
a=100*apply(a,1,sum)/genomesize
df=rbind(df, data.frame(percentage=a,type='PY'))
a=read.table("covered.0.05RPM.HC_nonNeuron.txt", header=F)$V2
a=100*sum(a)/genomesize
df=rbind(df, data.frame(percentage=a,type='NN'))
boxplot(percentage~type, df, outline =F, col=c('#F22A7B','#3182bd','#513931'), ylab="Covered percentage (%)", main="random 7 samples", yaxt="n")
axis(2, at=seq(15,30,5),labels=seq(15,30,5))
t.test(df$percentage[df$type=='SNDA'], mu=df$percentage[df$type=='NN'])
t.test(df$percentage[df$type=='PY'], mu=df$percentage[df$type=='NN'])
t.test(df$percentage[df$type=='SNDA'], df$percentage[df$type=='PY'])

# HCILB_SNDA with more details
df=read.table("covered.5reads.HCILB_SNDA.txt", header=F)
colnames(df)=c("samplecount", "coveredbp")
df=cbind(df, cumsum=cumsum(df$coveredbp))
par(mar=c(4,4,2,4), pty="s")
plot(df$cumsum, type='l', ylim=c(0,1900000000), xlim=c(1,90),ylab="", main="covered.5reads.HCILB_SNDA", yaxt="n", xaxt='n', xlab='')
points(df$coveredbp, type='h', lwd=4, col=colorRampPalette(c('blue','red'))(100), lend=2)
abline(h=df$coveredbp[1], lty=2)
axis(1, c(1,seq(10,90,10)), c(1,seq(10,90,10)))
axis(2, seq(0,2,.2)*1e9, labels=format(seq(0,2,.2),2), las=2)
axis(4, genomesize*seq(0,60,10)/100, labels=paste0(seq(0,60,10)), las=2)
mtext("Covered percentage (%)", 4, line=2)
mtext("Sample count", 1, line=2)
mtext("Covered base pairs (in billion)",2,line=3)

## aggregation plot for gene body coverage (to show degradation)
df=read.table("geneBodyCoverage.HCILB_SNDA.txt", header=F)
df=df[,-1]
df = (df-apply(df, 1, min))/(apply(df, 1, max)-apply(df,1,min))
df=data.frame(x=1:ncol(df),mean=colMeans(df), sd=apply(df, 2, sd))
library(ggplot2)
p=ggplot(data = df, aes(x = x, y = mean)) + theme_bw()
p=p+geom_line(aes(y = mean),size = 2) 
#p=p+geom_ribbon(aes(ymax = mean + sd, ymin = mean - sd), alpha = 0.5, fill = "grey70") # add ribbon for mean+/-sd
p
dev.off()

## ===============================================
# barplot of transcribed % comparison between GENCODE vs. BRAINCODE
## ===============================================

# GENCODE meta-exons (v19)
EXON=system("cut -f1-3 $GENOME/Annotation/Genes/exons.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T) # 122000567  (3.89%)  -- all exons
CDS=system("cut -f1-3 $GENOME/Annotation/Genes/cds.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T) # 34966072 (1.11%)  -- all CDS (protein-coding exon)
PC_EXON=system("fgrep protein_coding___protein_coding $GENOME/Annotation/Genes/exons.bed | cut -f1-3 | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T) # 75255917  (2.40%)  -- all protein-coding exons
# intergenic
intergenic=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed | awk '{s+=($3-$2)}END{print s}'",intern = T)  # 650101182
# EXONs
exons=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed -v | intersectBed -a - -b $GENOME/Annotation/Genes/exons.meta.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T)  # 101113507
# introns
introns=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed -v | intersectBed -a - -b $GENOME/Annotation/Genes/exons.meta.bed -v | awk '{s+=($3-$2)}END{print s}'",intern = T) # 1123233053


pdf("coverage.barplot.pdf", paper='us',width=4, height=4)
df=cbind(as.integer(c(EXON,0,0)),as.integer(c(exons,introns,intergenic)))
rownames(df) = c('exons', 'introns', 'intergenic')
colnames(df)=c("GENCODE", "BRAINCODE");
par(mar=c(4,4,2,4))
d=barplot(as.matrix(df), ylim=c(0,genomesize), col=c('#9ecae1', '#fc9272','#fec44f'), border =NA, axes=F, ylab="Human genome base pairs (in billion)")
text(x=d, y=apply(df,2,sum),pos=3, offset=.2, paste0(round(100*apply(df,2,sum)/genomesize,1),"%"), cex=2)
axis(2, at=c(0:3)*1e9, labels=0:3)
legend("topleft",col=c('#9ecae1', '#fc9272','#fec44f'), rownames(df), bty='n', pch=15)

# only show protein-coding exon [JUST FOR CLEMENS]
df=cbind(0,as.integer(c(exons,introns,intergenic)))
df=rbind(as.integer(c(CDS,0)),df)
rownames(df) = c('CDS','exons', 'introns', 'intergenic')
colnames(df)=c("GENCODE", "BRAINCODE");
par(mar=c(4,4,2,4))
d=barplot(as.matrix(df), ylim=c(0,genomesize), col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), border =NA, axes=F, ylab="Human genome base pairs (in billion)")
text(x=d, y=apply(df,2,sum),pos=3, offset=.2, paste0(round(100*apply(df,2,sum)/genomesize,1),"%"), cex=2)
axis(2, at=c(0:3)*1e9, labels=0:3)
legend("topleft",col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), rownames(df), bty='n', pch=15)

df=cbind(0,as.integer(c(exons,introns,intergenic)))
df=rbind(as.integer(c(PC_EXON,0)),df)
rownames(df) = c('protein_coding_genes.exons','exons', 'introns', 'intergenic')
colnames(df)=c("GENCODE", "BRAINCODE");
par(mar=c(4,4,2,4))
d=barplot(as.matrix(df), ylim=c(0,genomesize), col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), border =NA, axes=F, ylab="Human genome base pairs (in billion)")
text(x=d, y=apply(df,2,sum),pos=3, offset=.2, paste0(round(100*apply(df,2,sum)/genomesize,1),"%"), cex=2)
axis(2, at=c(0:3)*1e9, labels=0:3)
legend("topleft",col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), rownames(df), bty='n', pch=15)


dev.off()
