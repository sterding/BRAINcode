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
#bg=/data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bedGraph
bg=/data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bedGraph 

# use only the region above the basal level
basalLevel=`tail -n1 $bg | cut -f2 -d'=' | cut -f1`
awk -vmin=$basalLevel '$4>=min' $bg > /tmp/bg;
bedGraphToBigWig /tmp/bg $GENOME/Annotation/Genes/ChromInfo.txt ${bg/bedGraph/aboveBasal.bw}

ln -fs ${bg/bedGraph/aboveBasal.bw} RNAseq.aboveBasal.bigwig

bigWigAverageOverBed -minMax RNAseq.aboveBasal.bigwig $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 /tmp/RNAseq.aboveBasal.bigwig.exons.tab
bigWigAverageOverBed -minMax RNAseq.aboveBasal.bigwig $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.introns.bed12 /tmp/RNAseq.aboveBasal.bigwig.introns.tab

# bigWigAverageOverBed
   #name - name field from bed, which should be unique
   #size - size of bed (sum of exon sizes
   #covered - # bases within exons covered by bigWig
   #sum - sum of values over all bases covered
   #mean0 - average over bases with non-covered bases counting as zeroes
   #mean - average over just covered bases
   #min - minimal value of all bases
   #max - minimal value of all bases

join -1 1 -2 1 <(sort -k1,1 /tmp/RNAseq.aboveBasal.bigwig.exons.tab) <(sort -k1,1 /tmp/RNAseq.aboveBasal.bigwig.introns.tab) | awk '($3+$10)>0' > RNAseq.aboveBasal.bigwig.exons.intron.tab

awk '{print $2+$9, $0}' RNAseq.aboveBasal.bigwig.exons.intron.tab | sed 's/___/ /g' | cut -f1,2,3,6- -d' ' | sort -k3,3 -k1,1nr | awk '{if(id!=$3) {print; id=$3;}}END{if(id!=$3) print}' | head



# splicing_ratio = (exon-intro) / (exon + intron) in mean0 (average over bases with non-covered bases counting as zeroes)
# only keep genes with intron
join -1 1 -2 1 <(sort -k1,1 /tmp/RNAseq.aboveBasal.bigwig.exons.tab) <(sort -k1,1 /tmp/RNAseq.aboveBasal.bigwig.introns.tab) | awk '{OFS="\t"; splicing_ratio=(($5+$10)==0)?0:(($5-$10)/($5+$10)); rpkm_exon=$5*1000; rpkm_intron=$10*1000; cov_exon=$3/$2; cov_intron=$8/$7; print $1, rpkm_exon, rpkm_intron, cov_exon, cov_intron, splicing_ratio;}'>  /tmp/RNAseq.aboveBasal.bigwig.exons-introns.tab

echo -e "#chr\tstart\tend\tgeneID\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > gencode.v19.Tx.exons-introns.RNAseq.aboveBasal.bigwig.bed
cut -f1-4 $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 | sort -k4,4 | join -1 4 -2 1 - /tmp/RNAseq.aboveBasal.bigwig.exons-introns.tab -o '1.1,1.2,1.3,0,2.2,2.3,2.4,2.5,2.6' | sed 's/ /\t/g' >> gencode.v19.Tx.exons-introns.RNAseq.aboveBasal.bigwig.bed

# only the longest transcript per gene
echo -e "#chr\tstart\tend\tgeneID\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed
awk '{OFS="\t"; split($4,a,"___"); split(a[4],b,"."); $4=a[1]"___"a[2]"___"b[1]; if($1!~/^#/) print $0, $3-$2}' gencode.v19.Tx.exons-introns.RNAseq.aboveBasal.bigwig.bed | sort -k4,4 -k10,10nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-9 >> gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed

#awk '{if(($8/$7)>0.5 && $1 ~ /\.protein_coding/) print $0, $8/$7}' RNAseq.aboveBasal.bigwig.exons-introns.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt
#awk '{if(($8/$7)>0.5 && $10*10>=$5) print $0, $8/$7}' RNAseq.aboveBasal.bigwig.exons-introns.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt


############ R script ##########
R
df=read.table("gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed", header=F)
colnames(df)=c("chr","start","end","ID", "rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio") 

df0=subset(df, (rpkm_exon*rpkm_intron)>0 & grepl("protein_coding",ID))

attach(df0)



pdf("gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.pdf")
library(hexbin)
library(lattice)
x=log10(1e-6 + df0[,c('rpkm_exon', 'rpkm_intron')])
hexbinplot(rpkm_intron ~ rpkm_exon, data=x[rowMeans(x) > -6, ],
        xbins=80, asp=1,
        xlab="rpkm_exon in log10",
        ylab="rpkm_intron in log10",
        panel=function(...){
               panel.hexbinplot(...)
               panel.abline(a=0, b=1, col='black', lty=2, lwd=1)
        })

#my.colors <- function (n) { rev(heat.colors(n)) }
#hexbinplot(rpkm_intron ~ rpkm_exon,  xbins=80, data=x[rowMeans(x) > -6, ], asp=1, colramp = my.colors, colorcut = seq(0, 1, length = 20))

plot(x,
    pch=21,
    col=ifelse(cov_intron>0.3, '#ff0000',"#aaaaaa"),
    bg=ifelse(cov_exon>0.95, "#ff000033","#00000033"),
    xlab="rpkm_exon in log10", ylab="rpkm_intron in log10",
    xlim=range(-6.5, x), ylim=range(-6.5,x), xaxs = 'i', yaxs = 'i',
    cex=2*((end-start)-min(end-start))/(max((end-start))-min((end-start)))+0.3, asp=1);
abline(a=0, b=1, lty=2)
legend("topleft", c("intronic coverage > 30%", "intronic coverage <= 30%"), pch=21, col=c('#ff0000',"#aaaaaa"), pt.bg="#00000033", bty='n')

plot(cov_exon~cov_intron, df0, pch=20, cex=pmin(rpkm_intron,100)/100, col=rgb(0,0,1,pmin(rpkm_exon,800)/800))
abline(v=0.3, h=.4, lty=2, col='red', lwd=1)

hist(cov_intron[cov_exon>0.4], breaks=50)
hist(cov_exon[cov_intron>0.3], breaks=50)
        
hexbinplot(cov_intron ~ log10(end-start), asp=1)
dev.off()