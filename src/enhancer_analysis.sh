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
