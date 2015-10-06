## intron length vs. n_HITNE in introns
inputbed=~/eRNAseq/eRNA.bed
intersectBed -a $ANNOTATION/introns.meta.bed -b $inputbed -c | awk '{OFS="\t"; $5=($3-$2); print}' | sort -k5,5nr > introns.meta.nHTNE.txt

#R
df=read.delim("~/eRNAseq/introns.meta.nHTNE.txt", header=F)
colnames(df)=c('chr','start','end','id','length','strand','nHTNE')
head(df)
require(zoo)
plot(rollapply(df$length,width=200,by=40,FUN=mean),rollapply(df$nHTNE,width=200,by=40,FUN=mean))
