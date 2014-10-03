####################################
# Script to analyze active transcribed genes using RNAseq data
# Authors: Bioinformatics Team @ Scherzer's lab 
# Email: xdong@rics.bwh.harvard.edu
# Date: 9/4/2014
# Version: 1.0
####################################
#!/bin/bash

## FPKM for "gene" made of only introns
#######################################

cd ~/xd010/eRNAseq

# generate bigwig above the basal line
bg=/data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bedGraph
basalLevel=`tail -n1 $bg | cut -f2 -d'=' | cut -f1 -d' '`
awk -vmin=$basalLevel '$4>=min' $bg > /tmp/bg;
bedGraphToBigWig /tmp/bg $GENOME/Annotation/Genes/ChromInfo.txt ${bg/bedGraph/aboveBasal.bw}

ln -s /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.aboveBasal.bw RNAseq.aboveBasal.bigwig

bigWigAverageOverBed RNAseq.aboveBasal.bigwig $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 RNAseq.aboveBasal.bigwig.exons.tab
bigWigAverageOverBed RNAseq.aboveBasal.bigwig $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.introns.bed12 RNAseq.aboveBasal.bigwig.introns.tab

# splicing_ratio = exon / (exon + intron)
join -1 1 -2 1 <(sort -k1,1 RNAseq.aboveBasal.bigwig.exons.tab) <(sort -k1,1 RNAseq.aboveBasal.bigwig.introns.tab) | awk '{OFS="\t"; splicing_ratio=(($4+$9)==0)?0:($4/($4+$9)); rpkm_exon=$5*1000; rpkm_intron=$10*1000; cov_exon=$3/$2; cov_intron=$8/$7; print $1, rpkm_exon, rpkm_intron, cov_exon, cov_intron, splicing_ratio;}'>  RNAseq.aboveBasal.bigwig.exons-introns.tab  # only keep genes with intron

echo -e "#chr\tstart\tend\tgeneID\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > RNAseq.aboveBasal.bigwig.exons-introns.bed
cut -f1-4 $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 | sort -k4,4 | join -1 4 -2 1 - RNAseq.aboveBasal.bigwig.exons-introns.tab -o '1.1,1.2,1.3,0,2.2,2.3,2.4,2.5,2.6' |sed 's/ /\t/g' >> RNAseq.aboveBasal.bigwig.exons-introns.bed

#awk '{if(($8/$7)>0.5 && $1 ~ /\.protein_coding/) print $0, $8/$7}' RNAseq.aboveBasal.bigwig.exons-introns.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt
#awk '{if(($8/$7)>0.5 && $10*10>=$5) print $0, $8/$7}' RNAseq.aboveBasal.bigwig.exons-introns.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt


############ R script ##########
R
df=read.table("RNAseq.aboveBasal.bigwig.exons-introns.bed", header=F)
colnames(df)=c("chr","start","end","ID", "rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio") 

df0=df[grep("protein_coding.protein_coding",df$ID),]

attach(df0)



pdf("RNAseq.aboveBasal.bigwig.exons-introns.pdf")
library(hexbin)
x=log10(1e-6 + df0[,c('rpkm_exon', 'rpkm_intron')])
hexbinplot(rpkm_intron ~ rpkm_exon, data=x[rowMeans(x) > -6, ],
        xbins=80, asp=1,
        xlab="rpkm_exon in log10",
        ylab="rpkm_intron in log10",
        panel=function(x, y, ...){
               panel.hexbinplot(x,y,...)
               panel.abline(a=0, b=1, col='black', lty=2, lwd=1)
        })

#my.colors <- function (n) { rev(heat.colors(n)) }
#hexbinplot(rpkm_intron ~ rpkm_exon,  xbins=80, data=x[rowMeans(x) > -6, ], asp=1, colramp = my.colors, colorcut = seq(0, 1, length = 20))

plot(x,
    pch=21,
    col=ifelse(cov_intron>0.5, '#ff0000',"#aaaaaa"),
    bg=ifelse(cov_exon>0.95, "#ff000033","#00000033"),
    xlab="rpkm_exon in log10", ylab="rpkm_intron in log10",
    xlim=range(x), ylim=range(x),
    cex=(log(end-start)-min(log(end-start)))/(max(log(end-start))-min(log(end-start)))+0.0001, asp=1);
abline(a=0, b=1, lty=2)
legend("topleft", c("intronic coverage > 50%", "intronic coverage <= 50%"), pch=21, col=c('#ff0000',"#aaaaaa"), pt.bg="#00000033", bty='n')

hist(cov_intron[cov_exon>0.95], breaks=50)
hist(cov_exon[cov_intron>0.5], breaks=50)
        
hexbinplot(cov_intron ~ log10(end-start), asp=1)
dev.off()