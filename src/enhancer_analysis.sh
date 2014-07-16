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

# bigwig above the basal line
bg=/data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bedGraph
basalLevel=`tail -n1 $bg | cut -f2 -d'=' | cut -f1 -d' '`
awk -vmin=$basalLevel '$4>=min' $bg > /tmp/bg;
bedGraphToBigWig /tmp/bg $GENOME/Annotation/Genes/ChromInfo.txt ${bg/bedGraph/aboveBasal.bw}

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

ln -s /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.aboveBasal.bw RNAseq.aboveBasal.bigwig

bigWigAverageOverBed RNAseq.aboveBasal.bigwig /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v13.annotation.gtf.exons.bed12 gencode.v13.annotation.gtf.exons.bed12.RNAseq.bigwig.tab &
bigWigAverageOverBed RNAseq.aboveBasal.bigwig /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v13.annotation.gtf.introns.bed12 gencode.v13.annotation.gtf.introns.bed12.RNAseq.bigwig.tab &

join -1 1 -2 1 <(sort -k1,1 gencode.v13.annotation.gtf.exons.bed12.RNAseq.bigwig.tab) <(sort -k1,1 gencode.v13.annotation.gtf.introns.bed12.RNAseq.bigwig.tab) >  gencode.v13.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab

awk '{if(($8/$7)>0.5 && $1 ~ /\.protein_coding/) print $0, $8/$7}' gencode.v13.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt
awk '{if(($8/$7)>0.5 && $10*10>=$5) print $0, $8/$7}' gencode.v13.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab | tr -s "_" '\t'  | cut -f1-2 | uniq > genelist.cov0.5.pc.txt

############ R script ##########
R
df=read.table("gencode.v13.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.tab", header=F)
rownames(df)=df[,1];
df=df[,-1];
colnames(df)=c(paste("exons", c('size', 'covered', 'sum', 'mean0','mean'), sep="_"), paste("introns", c('size', 'covered', 'sum', 'mean0','mean'), sep="_")) 
attach(df)

pdf("gencode.v13.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.pdf")
library(hexbin)
x=log10(1e-6 + df[,c('exons_mean0', 'introns_mean0')]*1000)
hexbinplot(introns_mean0 ~ exons_mean0,  xbins=80, data=x[rowMeans(x) > -6, ],
        panel=function(x, y, ...){
               panel.hexbinplot(x,y,...)
               panel.abline(a=0, b=1, col='black', lty=2, lwd=1)
        })
xlab="exons_mean0 in log10(FPKM)",
        ylab="introns_mean0 in log10(FPKM)",
#my.colors <- function (n) { rev(heat.colors(n)) }
#hexbinplot(introns_mean0 ~ exons_mean0,  xbins=80, data=x[rowMeans(x) > -6, ], colramp = my.colors, colorcut = seq(0, 1, length = 20))

cov=introns_covered/introns_size
plot(x, pch=21, col=ifelse(cov>0.5, '#ff0000',"#aaaaaa"), bg="#00000033", xlab="exons_mean0 in log10(FPKM)", ylab="introns_mean0 in log10(FPKM)", xlim=range(x), ylim=range(x), cex=cov+0.0001, asp=1);
abline(a=0, b=1, lty=2)
legend("topleft", c("intronic coverage > 50%", "intronic coverage <= 50%"), pch=21, col=c('#ff0000',"#aaaaaa"), pt.bg="#00000033", bty='n')

dev.off()