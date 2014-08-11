## script to analyse eRNA defined by RNA-seq

cd ~/projects/PD/results/eRNA/externalData/RNAseq

# eRNA definition: 1) density 2x higher than the basal level, 2) summit >0.05 RPM and 3) located in non-generic regions (e.g. 500bp away from any exons)
for i in /data/neurogen/rnaseq_PD/results/merged/*.trimmedmean.uniq.normalized.bedGraph;
do
    basalLevel=`tail -n1 $i | cut -f2 -d'=' | cut -f1 -d' '`
    echo $i, $basalLevel;
    j=`basename ${i/bedGraph/eRNA.bed}`
    awk -vmin=$basalLevel '{OFS="\t"; if($4>=2*min) print $1,$2,$3,".",$4}' $i | mergeBed -scores max | awk '{OFS="\t"; if($4>=0.05) print $1,$2,$3,".",$4}' | mergeBed -d 100 -scores max | intersectBed -a - -b ../toExclude.bed -v > $j &
done

awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3}' HC_SNDA.trimmedmean.uniq.normalized.eRNA.bed > eRNA.bed

# measure the eRNA expression level (sum, mean0, intron/exon rate of host gene (if any))
for i in ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw; do
    echo $i;
    bigWigAverageOverBed $i eRNA.bed stdout | cut -f1,5 | sort -k1,1 > $i.eRNA.measured 
done

# sed 's/_/\t/g' | cut -f1-3,7 | intersectBed -a stdin -b gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.bed -wao | awk '{OFS="\t"; print $0, ($15==".")?-1:(0+$15/$14); }' | groupBy -g 1,2,3,4 -c 20 -o max 


# Version2: eRNA definition: 1) p-value <0.05 more enriched than the background (non-genic regions), 2) summit >0.05 RPM and 3) located in non-generic regions (e.g. 500bp away from any exons)
#1: create random regions (100,000 regions of 300bp each) for background
#1a: blacklist region to exclude from background
ANNOTATION=$GENOME/Annotation/Genes
cat $ANNOTATION/gencode.v19.annotation.bed12 $ANNOTATION/knownGene.bed12 $ANNOTATION/NONCODEv4u1_human_lncRNA.bed12 | bed12ToBed6 | cut -f1-3 | grep -v "_" |slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 > /tmp/bg.bed
cut -f1-3 $ANNOTATION/rRNA.bed >> /tmp/bg.bed  # rRNA
grep -v track ~/projects/PD/results/eRNA/externalData/CAGE/permissive_enhancers.bed| cut -f1-3 >> /tmp/bg.bed # CAGE-enhancer
#cat $ANNOTATION/SINE.bed $ANNOTATION/LINE.bed | cut -f1-3 >> /tmp/bg.bed  # LINE and SINE
cat /tmp/bg.bed | sortBed | mergeBed > blacklist.bed

#1b: 
grep -v chrM $GENOME/Annotation/Genes/ChromInfo.txt > ChromInfo.nochrM.txt
ls ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw ~/neurogen/rnaseq_PD/results/merged/*trimmedmean.uniq.normalized.bw | \
    parallel 'echo {};  bedtools random -g ChromInfo.nochrM.txt -l 300 -n 500000 | bedtools intersect -a - -b blacklist.bed -v | head -n100000 | bigWigAverageOverBed {} stdin {}.rdbg'

#2: distribution of random background, in order to define the cutoff with p=0.0001 significance

R
path="~/neurogen/rnaseq_PD/results/merged/*trimmedmean.uniq.normalized.bw.rdbg"; filename="background_merged.pdf"

path="~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw.rdbg"; filename="background.pdf"

pdf(filename)
DF=c()
Fn0;
for(i in Sys.glob(path)){
    print(i)
    df=read.table(i, header=F)[,5]
    df[df>2]=2
    DF=c(DF, df)
    Fn=ecdf(sample(df, 100000))
    if(grepl("HC_SNDA.trimmedmean", i)) Fn0=Fn;
    main=ifelse(grepl("merged", i), sub(".*merged/(.*).bw.*","\\1", i), sub(".*run_output/(.*)/uniq.*","\\1", i))
    plot(Fn, verticals = TRUE, do.points = FALSE, main=main, ylim=c(0.99, 1), xlab="mean RPM", ylab="cummulative percentage (approx. 1-p)")
    inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
    abline(h=0.999, v=g(0.999), col='red', lty=2, lwd=1)
    points(g(0.999), 0.999, col='red', pch=19)
    #text(g(0.999), 0.999, paste(" (",round(g(0.999),2), ", 99.9%)", sep=""), cex=3, adj=c(0,1))
    text(g(0.999), 0.999, round(g(0.999),2), cex=5, adj=c(0,1))
}
# all samples
Fn=ecdf(sample(DF, 100000))
plot(Fn, verticals = TRUE, do.points = FALSE, main="All samples", ylim=c(0.99, 1), xlab="mean RPM", ylab="cummulative percentage (approx. 1-p)")
inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
abline(h=0.999, v=g(0.999), col='red', lty=2, lwd=1)
points(g(0.999), 0.999, col='red', pch=19)
text(g(0.999), 0.999, round(g(0.999),2), cex=5, adj=c(0,1))

dev.off()

# significance
path="~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw"
EXP=data.frame()
PV=data.frame()
QV=data.frame()
for(i in Sys.glob(path)){
    print(i)
    
    # read background
    df=read.table(paste(i,"rdbg",sep="."), header=F)[,5]
    Fn=ecdf(df)
    
    # read expression
    expression=read.table(paste(i,"eRNA.measured",sep="."), header=F)
    pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));
    qvalue=as.numeric(format(p.adjust(pvalue, "BH"), digits=3));
    write.table(cbind(expression[,1:2], pvalue=pvalue, qvalue=qvalue), file=paste(i,"eRNA.measured.significance",sep="."), quote=F, sep ="\t", col.names =F, row.names=F)
    
    # merge
    if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; expression[,2]=qvalue; QV=expression; }
    else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); QV=cbind(QV, qvalue); }
}
        


# add p-values for eRNA regions
eRNA=read.table("HC_SNDA.trimmedmean.uniq.normalized.eRNA.measured", header=F)
colnames(eRNA)=c('chr','start', 'end', 'summit', 'sum', 'mean')
eRNA=cbind(eRNA, pvalue=1-Fn0(eRNA$mean), qvalue=round(p.adjust(1-Fn0(eRNA$mean), "BH"), 3))
write.table(eRNA, "HC_SNDA.trimmedmean.uniq.normalized.eRNA.measured.xls", sep="\t", quote = F, col.names = T, row.names = F)

pdf("HC_SNDA.trimmedmean.uniq.normalized.eRNA.measured.pdf")
hist(eRNA$pvalue, breaks=50, main="p-value")
hist(eRNA$qvalue, breaks=50, main="BenjaminiÐ Hochberg adjusted p-value")
dev.off()

# RPM with p=0.001 in all samples is: 0.64
# RPM with p=0.001 in HC_SNDA.trimmedmean.uniq.normalized is: 0.17

#3: robust set of eRNA


# length distribution
for i in *.trimmedmean.uniq.normalized.eRNA.bed; do echo $i; wc -l $i; awk '{print $3-$2}' $i | textHistogram -binSize=20 -maxBinCount=50 stdin; done

awk '{OFS="\t"; if(($3-$2)>=200) print $1, $2, $3, $1"_"$2"_"$3, $4}' HC_SNDA.trimmedmean.uniq.normalized.eRNA.bed > eRNA.bed

# Total counts of CAGE reads
toBinRegionsOnBigwig.sh ../CAGE/ctssTotalCounts.fwd.bw eRNA.bed 1 max > eRNA.CAGE.fwd.bed &
toBinRegionsOnBigwig.sh ../CAGE/ctssTotalCounts.rev.bw eRNA.bed 1 max > eRNA.CAGE.rev.bed &

# TF count
toBinRegionsOnBigwig.sh ../TFBS/TFBS.bigwig eRNA.bed 1 max > eRNA.TFBS.bed &

# Histone
toBinRegionsOnBigwig.sh ../Histone/Histone.SN.H3K27ac.bigwig eRNA.bed 1 max > eRNA.SN.H3K27ac.bed &

# DNase
toBinRegionsOnBigwig.sh ../DNase/DNase.bigwig eRNA.bed 1 max > eRNA.DNase.bed &

echo -e "position\tRNAseq\tCAGE.fwd\tCAGE.rev\tDNase\tH3K27ac\tTFBS" > eRNA_merged.txt
paste eRNA.*bed | sed 's/ /\t/g' | cut -f4,5,7,9,11,13,15 >> eRNA_merged.txt

R
df=read.table("eRNA_merged.txt", header=T)
rownames(df)=df[,1]; df=df[,-1]
attach(df)
x=log10(1+cbind(CAGE.fwd, CAGE.rev))
plot(x, pch=1, col='#ff000022', cex=sqrt(RNAseq))

library(flashClust)
d=df[,grep('CAGE',colnames(df))]
d=d[sample(nrow(d),2000),]
dis=dist(d)
hc <- hclust(dis)
plot(hc)

df=df[with(df, order(RNAseq)),]
for(i in 1:ncol(df)){
    image(t(df[,i,drop=F]))
}


base_dir=/apps/source/mpich2_1.4.1
source ${base_dir}/mpich2lsf.sh $LSB_HOSTS

mpirun  -np $nproc  -machinefile ${LSB_JOBID}_machines


## FPKM for "gene" made of only introns
#######################################
# bigwig above the basal line
bg=/data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bedGraph
basalLevel=`tail -n1 $bg | cut -f2 -d'=' | cut -f1 -d' '`
awk -vmin=$basalLevel '$4>=min' $bg > /tmp/bg;
bedGraphToBigWig /tmp/bg $GENOME/Annotation/Genes/ChromInfo.txt ${bg/bedGraph/aboveBasal.bw}

ln -s /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.aboveBasal.bw RNAseq.aboveBasal.bigwig

bigWigAverageOverBed RNAseq.aboveBasal.bigwig /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 gencode.v19.annotation.gtf.exons.bed12.RNAseq.bigwig.tab &
bigWigAverageOverBed RNAseq.aboveBasal.bigwig /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.gtf.introns.bed12 gencode.v19.annotation.gtf.introns.bed12.RNAseq.bigwig.tab &

join -1 1 -2 1 <(sort -k1,1 gencode.v19.annotation.gtf.exons.bed12.RNAseq.bigwig.tab) <(sort -k1,1 gencode.v19.annotation.gtf.introns.bed12.RNAseq.bigwig.tab) >  gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab

cut -f1-4 /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.gtf.exons.bed12 | sort -k4,4 | join -1 4 -2 1 - gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab -o '1.1,1.2,1.3,0,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11' |sed 's/ /\t/g' > gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.bed

#awk '{if(($8/$7)>0.5 && $1 ~ /\.protein_coding/) print $0, $8/$7}' gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt
#awk '{if(($8/$7)>0.5 && $10*10>=$5) print $0, $8/$7}' gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt



############ R script ##########
R
df=read.table("gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab", header=F)
rownames(df)=df[,1];
df=df[,-1];
colnames(df)=c(paste("exons", c('size', 'covered', 'sum', 'mean0','mean'), sep="_"), paste("introns", c('size', 'covered', 'sum', 'mean0','mean'), sep="_")) 
attach(df)

pdf("gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.pdf")
library(hexbin)
x=log10(1e-6 + df[,c('exons_mean0', 'introns_mean0')]*1000)
hexbinplot(introns_mean0 ~ exons_mean0, data=x[rowMeans(x) > -6, ],
        xbins=80, asp=1,
        xlab="exons_mean0 in log10(FPKM)",
        ylab="introns_mean0 in log10(FPKM)",
        panel=function(x, y, ...){
               panel.hexbinplot(x,y,...)
               panel.abline(a=0, b=1, col='black', lty=2, lwd=1)
        })

#my.colors <- function (n) { rev(heat.colors(n)) }
#hexbinplot(introns_mean0 ~ exons_mean0,  xbins=80, data=x[rowMeans(x) > -6, ], asp=1, colramp = my.colors, colorcut = seq(0, 1, length = 20))

cov=introns_covered/introns_size
plot(x, pch=21, col=ifelse(cov>0.5, '#ff0000',"#aaaaaa"), bg="#00000033", xlab="exons_mean0 in log10(FPKM)", ylab="introns_mean0 in log10(FPKM)", xlim=range(x), ylim=range(x), cex=cov+0.0001, asp=1);
abline(a=0, b=1, lty=2)
legend("topleft", c("intronic coverage > 50%", "intronic coverage <= 50%"), pch=21, col=c('#ff0000',"#aaaaaa"), pt.bg="#00000033", bty='n')

dev.off()

eRNA=read.table("HC_SNDA.trimmedmean.uniq.normalized.eRNA.measured.introCov", header=F)
colnames(eRNA)=c('chr','start','end', 'summit', 'size', 'covered', 'sum', 'mean0','mean','hostgene_intron_cov')
