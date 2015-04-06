## script to analyse eRNA defined by HCILB RNA-seq

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

cd ~/projects/PD/results/eRNA/externalData/RNAseq

inputBG=/data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bedGraph

################################################
# eRNA definition:
# 1) density higher than the basal level,  
# 2) summit >0.05 RPM, --> p<0.05 comparing to the transcriptional noise
# 3) located in non-generic regions (e.g. 500bp away from any annotated exons),
# 4) at least 100bp in length,
# 5) not from highly actively transcribing genes (e.g. pre-mRNA, intronic coverage > 50% & exonic coverage >90%)
# 6) q-value<0.05 when comparing with random non-functional background
################################################

#: background region to measure transcriptional noise: genomic regions excluding the known regions with RNA activities (exons+/-500bp, rRNA, CAGE-defined enhancers, promoters)
ANNOTATION=$GENOME/Annotation/Genes
cat $ANNOTATION/gencode.v19.annotation.bed12 $ANNOTATION/knownGene.bed12 $ANNOTATION/NONCODEv4u1_human_lncRNA.bed12 | bed12ToBed6 | cut -f1-3 | grep -v "_" | slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 > /tmp/bg.bed
cut -f1-3 $ANNOTATION/rRNA.bed >> /tmp/bg.bed  # rRNA
grep -v track ~/projects/PD/results/eRNA/externalData/CAGE/TSS_human.bed | grep -v "211,211,211" | cut -f1-3 | grep -v "_" | slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 >> /tmp/bg.bed # +/-500bp flanking around the CAGE-predicted TSS (downloaded from: http://fantom.gsc.riken.jp/5/datafiles/latest/extra/TSS_classifier/)
#cat $ANNOTATION/SINE.bed $ANNOTATION/LINE.bed | cut -f1-3 >> /tmp/bg.bed  # LINE and SINE
cat /tmp/bg.bed | sortBed | mergeBed > ../toExclude.bed
grep -v track ~/projects/PD/results/eRNA/externalData/CAGE/permissive_enhancers.bed | cut -f1-3 >> /tmp/bg.bed # CAGE-enhancer
cat /tmp/bg.bed | sortBed | mergeBed > ../blacklist.bed

# RNAseq signal distribution in the background region
intersectBed -a $inputBG -b ../blacklist.bed -sorted -v | awk '{OFS="\t"; print $3-$2, $4}' | shuf > transcriptional.noise.rpm.txt

#R
df=read.table("transcriptional.noise.rpm.txt", comment.char = "", nrows = 2000000)
df=log10(as.numeric(do.call('c',apply(df, 1, function(x) rep(x[2],x[1])))))
library(fitdistrplus)
fitn=fitdist(df,'norm')
pdf("transcriptional.noise.distribution.pdf", width=8, height=6)
hist(df, breaks=100, prob=TRUE, xlab='log10(RPM)', main='Distribution of transcriptional noise')
lines(density(df, bw=0.15))
m=round(as.numeric(fitn$estimate[1]),digits=3)
sd=round(as.numeric(fitn$estimate[2]),digits=3)
lines(density(rnorm(n=2000000, mean=m, sd=sd),bw=0.25), col='blue',lty=2)
p=round(qnorm(.05, mean=m, sd=sd, lower.tail = F), digits=3)
lines(y=c(0,0.3),x=c(p,p),col='red')
text(p,0.2,paste0("P(X>",p,") = 0.05\nRPM=10**",p,"=",round(10**p,digits=3)), adj=c(0,0))
legend("topright", c("empirical density curve", paste0("fitted normal distribution \n(mean=",m,", sd=",sd,")")), col=c('black','blue'), lty=c(1,2), bty='n')
dev.off()

# Dsig: 10**-1.105 == 0.079

## any region with RPM density > 0.101
#basalLevel=0.101
#j=`basename ${inputBG/bedGraph/eRNA.bed}`
#awk -vmin=$basalLevel '{OFS="\t"; if($4>min) print $1,$2,$3,".",$4}' $inputBG | mergeBed -d 100 -scores max | intersectBed -a - -b ../toExclude.bed -v > $j
##wc -l $j
##40451 trimmedmean.uniq.normalized.HCILB_SNDA.eRNA.bed

#for i in /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bedGraph;
#do
#    basalLevel=`tail -n1 $i | cut -f2 -d'=' | cut -f1 -d' '`
#    echo $i, $basalLevel;
#    j=`basename ${i/bedGraph/eRNA.bed}`
#    awk -vmin=$basalLevel '{OFS="\t"; if($4>=2*min) print $1,$2,$3,".",$4}' $i | mergeBed -scores max | awk '{OFS="\t"; if($4>=0.05) print $1,$2,$3,".",$4}' | mergeBed -d 100 -scores max | intersectBed -a - -b ../toExclude.bed -v > $j &
#done

# step1: any regions with summit RPM > peakLevel and border > baseLevel
# =================
basalLevel=`tail -n1 $inputBG | cut -f2 -d'=' | cut -f1`
awk -vmin=$basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' $inputBG | mergeBed -scores max > eRNA.tmp1

# step2: summit RPM >=Dsig (density with p<0.05)
# =================
Dsig=0.079
awk -vD=$Dsig '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' eRNA.tmp1 | mergeBed -d 100 -scores max > eRNA.tmp2

# step3: located in non-generic regions (e.g. 500bp away from any annotated exons),
# =================
intersectBed -a eRNA.tmp2 -b ../toExclude.bed -v > eRNA.tmp3

# step4: length > 100nt
# =================
awk '{OFS="\t"; if(($3-$2)>100) print $1,$2,$3,$1"_"$2"_"$3}' eRNA.tmp3 > eRNA.tmp4

# step5: calculate the significance of eRNA
# =================
#1: create 100,000 random regions (400bp each) as background and calculate their signals
for i in ~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.normalized.bw ~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bw;
do
    bsub -q normal -n 1 "bedtools random -g $GENOME/Annotation/Genes/ChromInfo.txt -l 400 -n 500000 | grep -v chrM | intersectBed -a - -b ../blacklist.bed -v | head -n100000 | bigWigAverageOverBed $i stdin $i.rdbg";
    bsub -q normal -n 1 "bigWigAverageOverBed $i eRNA.tmp4 stdout | cut -f1,5 | sort -k1,1 | awk '{OFS=\"\t\"; print \$1, \$2*1000+0}'> $i.eRNA.RPKM"
done

### 2: distribution of random background, in order to define the cutoff with p=0.0001 significance
R
# significance
path=c("~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bw", "~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.normalized.bw")
pdf("background.RNAseq.cummulative.plot.pdf")

## read in only the 80 subjects/samples (w/ genotype)
IDs=read.table('~/neurogen/rnaseq_PD/results/merged/RNAseqID.wGenotyped.list',stringsAsFactors =F)[,1]
IDs=IDs[grep("^[HI].*_SNDA", IDs)]

EXP=data.frame(); PV=data.frame(); QV=data.frame(); id="locus"
for(i in Sys.glob(path)){
    ii=ifelse(grepl("merged", i), sub(".*merged/(.*).bw.*","\\1", i), sub(".*run_output/(.*)/uniq.*","\\1", i));
    if(! (ii %in% IDs || grepl("trimmedmean",ii))) next;
    print(i)
    # read background
    df=read.table(paste(i,"rdbg",sep="."), header=F)[,5] * 1000  # convert mean RPM to RPKM
    Fn=ecdf(df)
    
    # plot the cummulative plot
    plot(Fn, verticals = TRUE, do.points = FALSE, main=ii, ylim=c(0.99, 1), xlab="RPKM", ylab="cummulative percentage (approx. 1-p)")
    inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
    abline(h=0.999, v=g(0.999), col='red', lty=2, lwd=1)
    points(g(0.999), 0.999, col='red', pch=19)
    text(g(0.999), 0.999, round(g(0.999),2), cex=5, adj=c(0,1))
    
    if(grepl("trimmedmean",ii)) next;
    id=c(id, ii)

    # read expression
    expression=read.table(paste(i,"eRNA.RPKM",sep="."), header=F)
    pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));
    qvalue=as.numeric(format(p.adjust(pvalue, "BH"), digits=3));
    write.table(cbind(expression[,1:2], pvalue=pvalue, qvalue=qvalue), file=paste(i,"eRNA.RPKM.significance",sep="."), quote=F, sep ="\t", col.names =F, row.names=F)
    
    # merge
    if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; expression[,2]=qvalue; QV=expression; }
    else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); QV=cbind(QV, qvalue); }
}
dev.off()

colnames(EXP)=id; colnames(PV)=id; colnames(QV)=id;
rM=rowMeans(QV[,-1]<=0.05)
write.table(EXP[rM>0.25,], "eRNA.80samples.RPKM.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(PV[rM>0.25,], "eRNA.80samples.pvalue.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(QV[rM>0.25,], "eRNA.80samples.qvalue.xls", col.names=T, row.names=F, sep="\t", quote=F)

pdf("eRNA.80samples.qvalue.hist.pdf", width=8, height=6)
h=hist(rM, breaks=80, xlim=c(0,1), main="",xlab=expression("Percentage of HC/ILB SNDA samples (out of 80) with q-value" <= "0.05"), ylab="Count of eRNAs", freq=T)
abline(v=0.25, lty=2, col='red')
legend(0.25, max(h$counts), c(bquote(.(sum(rM>0.25)) ~ "eRNAs"), expression("with q-value" <= "0.05"), "in at least 25% of samples"), adj=c(0,1), bty='n', text.col='red', cex=1.5)
dev.off()

q('no')
## R end

# final set of eRNAs
# =================
awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) print a[1],a[2],a[3],$1}' eRNA.80samples.RPKM.xls > eRNA.bed
rsync -azv eRNA.bed xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/merged

################################################################################################
# characterize eRNA
################################################################################################


#############################################################
# clustering of eRNA expression
#############################################################
Rscript eRNA.clustering.R eRNA.80samples.RPKM.xls
Rscript eRNA.clustering.R eRNA.80samples.rawcount.xls


#############################################################
# DESeq of eRNA (by BZ)
#############################################################



#############################################################
# target genes of eRNA assigned by correlation of expression
#############################################################

cut -f1,9- ~/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls > genes.fpkm.allSamples.uniq.xls
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R eRNA.allsamples.RPKM.tab genes.fpkm.allSamples.uniq.xls eRNA.correlate.gene.in.RPKM
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R eRNA.allsamples.readscount.tab ~/neurogen/rnaseq_PD/results/DE_DESeq2/PDvsHC/htseqcount.raw.allsamples.xls eRNA.correlate.gene.in.readscount.tab

# join pcc and rho
join -j 1 <(awk '{print $1"__"$2,$3;}' eRNA.correlate.gene.in.RPKM.pcc.tab|sort) <(awk '{print $1"__"$2,$3;}' eRNA.correlate.gene.in.RPKM.rho.tab|sort) | sort -r | sed 's/__/\t/g;s/ /\t/g' >eRNA.correlate.gene.in.RPKM.cor.tab

# eRNA per gene
cut -f2 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | datamash mean 1
#9.451491660794
# gene per eRNA
cut -f1 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | datamash mean 1
#15.613116026387

# top eRNA with most target genes
cut -f1 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sort -k1,1nr

# xyplot

grep chr eRNA.correlate.gene.in.RPKM.rho.tab | head -n20 | while read x y rest ; do
echo $x, $y;
set -v
    grep -P "$x|locus" eRNA.allsamples.RPKM.tab > /tmp/xyplot.tab
    grep $y genes.fpkm.allSamples.uniq.xls >> /tmp/xyplot.tab
    Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.xyplot.R /tmp/xyplot.tab "$x.$y"
set +v
done

#############################################################
# overlap of GWAS SNPs with eRNA, exons etc. 
#############################################################
awk 'BEGIN{FS="\t";} $9<=1e-8' /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/gwascatalog2014AUGID.selected.txt | cut -f2 | sed 's/ (.*//g' | sort | uniq -c | sort -k2,2 | awk '$1>30'




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
