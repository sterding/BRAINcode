## intron length vs. n_HITNE in introns
inputbed=~/eRNAseq/eRNA.bed
intersectBed -a $ANNOTATION/introns.meta.bed -b $inputbed -c | awk '{OFS="\t"; $5=($3-$2); print}' | sort -k5,5nr > introns.meta.nHTNE.txt
bedtools shuffle -excl $ANNOTATION/exons.meta.bed -noOverlapping -i $inputbed -g $ANNOTATION/ChromInfo.txt | intersectBed -a $ANNOTATION/introns.meta.bed -b stdin -c | awk '{OFS="\t"; $5=($3-$2); print}' | sort -k5,5nr > introns.meta.nHTNEshuffled.txt

## gene length vs. n_HITNE in genes
inputbed=~/eRNAseq/eRNA.bed
awk '{OFS="\t"; print $1,$2,$3,$4"___"$5"___"$6,0,$7}' $ANNOTATION/genes.bed | intersectBed -a stdin -b $inputbed -c | awk '{OFS="\t"; $5=($3-$2); print}' | sort -k5,5nr > genes.vs.nHTNE.txt
bedtools shuffle -excl $ANNOTATION/exons.meta.bed -noOverlapping -i $inputbed -g $ANNOTATION/ChromInfo.txt | intersectBed -a <(awk '{OFS="\t"; print $1,$2,$3,$4"___"$5"___"$6,0,$7}' $ANNOTATION/genes.bed ) -b stdin -c | awk '{OFS="\t"; $5=($3-$2); print}' | sort -k5,5nr > genes.vs.nHTNEshuffled.txt


#R
df=read.delim("~/eRNAseq/introns.meta.nHTNE.txt", header=F)
colnames(df)=c('chr','start','end','id','length','strand','nHTNE')
head(df)
require(zoo);library(ggplot2)

q <- data.frame(x=rollapply(df$length,width=200,by=40,FUN=mean),
                mean=rollapply(df$nHTNE,width=200,by=40,FUN=mean),
                sem=rollapply(df$nHTNE,width=200,by=40,FUN=function(z) sd(z)/sqrt(200)),
                type='HTNE')
df=read.delim("~/eRNAseq/introns.meta.nHTNEshuffled.txt", header=F)
colnames(df)=c('chr','start','end','id','length','strand','nHTNE')
q <- rbind(q, data.frame(x=rollapply(df$length,width=200,by=40,FUN=mean),
                mean=rollapply(df$nHTNE,width=200,by=40,FUN=mean),
                sem=rollapply(df$nHTNE,width=200,by=40,FUN=function(z) sd(z)/sqrt(200)),
                type='random'))
ggplot(data = q, aes(x = x/1000, y = mean, group=type)) + xlab("Intron length (kb)") + ylab("Mean HTNEs number") +
  geom_line(size=.8, aes(colour=factor(type))) + scale_color_manual(values=c("#ff0000", "#000000")) +
  geom_ribbon(aes(ymax = mean + sem, ymin = mean - sem, fill=factor(type)), alpha = 0.2) + scale_fill_manual(values=c("#ff0000", "#000000"))
ggsave("~/eRNAseq/introns.vs.nHTNE.pdf")


df=read.delim("~/eRNAseq/genes.vs.nHTNE.txt", header=F)
colnames(df)=c('chr','start','end','id','length','strand','nHTNE')
head(df)
require(zoo);library(ggplot2)

q <- data.frame(x=rollapply(df$length,width=200,by=40,FUN=mean),
                mean=rollapply(df$nHTNE,width=200,by=40,FUN=mean),
                sem=rollapply(df$nHTNE,width=200,by=40,FUN=function(z) sd(z)/sqrt(200)),
                type='HTNE')
df=read.delim("~/eRNAseq/genes.vs.nHTNEshuffled.txt", header=F)
colnames(df)=c('chr','start','end','id','length','strand','nHTNE')
q <- rbind(q, data.frame(x=rollapply(df$length,width=200,by=40,FUN=mean),
                         mean=rollapply(df$nHTNE,width=200,by=40,FUN=mean),
                         sem=rollapply(df$nHTNE,width=200,by=40,FUN=function(z) sd(z)/sqrt(200)),
                         type='random'))
ggplot(data = q, aes(x = x/1000, y = mean, group=type)) + xlab("Gene length (kb)") + ylab("Mean HTNEs number") +
  geom_line(size=.8, aes(colour=factor(type))) + scale_color_manual(values=c("#ff0000", "#000000")) +
  geom_ribbon(aes(ymax = mean + sem, ymin = mean - sem, fill=factor(type)), alpha = 0.2) + scale_fill_manual(values=c("#ff0000", "#000000"))

ggsave("~/eRNAseq/genes.vs.nHTNE.pdf")


