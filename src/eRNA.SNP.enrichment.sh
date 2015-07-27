# The Fisher exact test was performed to test the odds ratio of trait association with enhancers vs non-enhancers taking into account the genomic distribution of SNPs. 
# 
# (A/B) / (C/D), where
# 
# A: number of trait-associated SNPs in enhancers
# B: number of non-trait-associated SNPs in enhancers (from dbSNP - GWAS SNPs for that trait)
# C: number of trait-associated SNPs distal to enhancers
# D: number of non-trait-associated SNPs distal to enhancers (from dbSNP - GWAS SNPs for that trait)

# convert GWAS SNPs into diseases:SNPs_in_LD

#cat gwascatalog2015Apr.gwas-clean-hg19.uniq.bed | while IFS=$'\t' read chr start end name score strand disease rest; do echo -e $chr"\t"$start"\t"$end"\t"$name"\t"$score"\t"$strand"\t"$disease"\t"`fgrep -w $end ../Variation/1000G/LD/LD_$chr.ld | awk -vend=$end '{if($7>0.4) {id=($2==end)?$6:$3; ID=id";"ID; e=($2==end)?$5:$2; E=e";"E; R2=$7";"R2; D=$8";"D;;}}END{printf ID"|"E"|"R2"|"D}'`; done > gwascatalog2015Apr.gwas-clean-hg19.uniq.snps_in_LD.bed

snps_in_LD=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.snps_in_LD.r2.0.4.bed

# all
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.all
# eRNA 
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b eRNA.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.eRNA
# mRNA exons 
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq| intersectBed -a stdin -b <(grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed)  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.exon
# promoter -- [-300,+100] of TSS peak summit, as http://promoter.binf.ku.dk/documentation/ described
# curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz | zcat | awk '{OFS="\t"; s=($6=="+")?($7-300):($7-100); if(s<0) s=0; print $1,s,s+400}' | sortBed | mergeBed -i stdin | awk '{s+=($3-$2)}END{print s}'
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b <(curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz | zcat | awk '{OFS="\t"; s=($6=="+")?($7-300):($7-100); if(s<0) s=0; print $1,s,s+400}')  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.promoter

# randomly sampling
bedtools shuffle -excl blacklist.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > eRNA.random.bed
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b eRNA.random.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.random

# total SNPs in dbSNP v137: 
wc -l /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed.groupped.SNP
# 48709140
intersectBed -a /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.bed -u | wc -l
# 830459
grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed | intersectBed -a /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed.groupped.SNP -b stdin -u | wc -l
# 1767413
curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz | zcat | awk '{OFS="\t"; s=($6=="+")?($7-300):($7-100); if(s<0) s=0; print $1,s,s+400}' | intersectBed -a /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed.groupped.SNP -b stdin -u | wc -l
# 835088
intersectBed -a /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.random.bed -u | wc -l
# 861067

# R test
setwd("~/eRNAseq")
n1=830459; n2=48709140;
all=read.table("SNP.count.all"); rownames(all)=all[,1]
x=read.table("SNP.count.eRNA"); rownames(x)=x[,1]
df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
results=cbind(Disease_or_Trait=rownames(df), df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type='HiTNEs')

n1=1767413; n2=48709140;
x=read.table("SNP.count.exon"); rownames(x)=x[,1]
df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type='exons'))

n1=835088; n2=48709140;
x=read.table("SNP.count.promoter"); rownames(x)=x[,1]
df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type='promoter'))

n1=861067; n2=48709140;
x=read.table("SNP.count.random"); rownames(x)=x[,1]
df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type='random'))

results=subset(results, pvalue<0.001 & observed>3)
results = results[with(results, order(type, -OR)), ]
table(results$type)
# HiTNEs    exons promoter   random 
# 28       46       61       27 
results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
write.table(results, "eRNA.SNP.enrichments.xls", sep="\t", col.names = T, row.names = F)

require(ggplot2)
df=subset(results, type=='HiTNEs')
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, as.character(results$Disease_or_Trait))
ggplot(results, aes(x=Disease_or_Trait2, y=OR, fill=type)) + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=9), legend.justification=c(1,1), legend.position=c(1,1)) + geom_text(aes(label=paste0(observed," (",round(-log10(pvalue),1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=3) + ylim(0, 170)

ggsave("eRNA.SNP.enrichment.pdf")


## BELOW ARE CODE TO TEST per-bp enrichment test (WHICH IS NOT PROPER FOR THIS CASE)
# =========================================================================

# # all (3e9)
# awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.all
# # eRNA (51,446,599)
# awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b eRNA.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.eRNA
# # mRNA exons (75,255,917)
# awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq| intersectBed -a stdin -b <(grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed)  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.exon
# # promoter (45,683,963) -- [-300,+100] of TSS peak summit, as http://promoter.binf.ku.dk/documentation/ described
# # curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz | zcat | awk '{OFS="\t"; s=($6=="+")?($7-300):($7-100); if(s<0) s=0; print $1,s,s+400}' | sortBed | mergeBed -i stdin | awk '{s+=($3-$2)}END{print s}'
# awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b <(curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz | zcat | awk '{OFS="\t"; s=($6=="+")?($7-300):($7-100); if(s<0) s=0; print $1,s,s+400}')  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.promoter
# # random (51,446,599)
# ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
# sortBed -i $ANNOTATION/exons.meta.bed | mergeBed -i - | cat - $ANNOTATION/hg19.gap.bed  > /tmp/toexclude.bed
# # randomly sampling
# tmp=`mktemp`
# bedtools shuffle -excl /tmp/toexclude.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp
# awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD |  sortBed | uniq | intersectBed -a stdin -b $tmp -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.count.random
# 
# # R test
# setwd("~/eRNAseq")
# all=read.table("SNP.count.all"); rownames(all)=all[,1]
# x=read.table("SNP.count.eRNA"); rownames(x)=x[,1]
# df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
# results=cbind(Disease_or_Trait=rownames(df), df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],51,x[2]-x[1], 3137-51), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],51,x[2]-x[1], 3137-51), nrow = 2), alternative='greater')$estimate), type='HiTNEs')
# 
# x=read.table("SNP.count.exon"); rownames(x)=x[,1]
# df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
# results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],75,x[2]-x[1], 3137-75), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],75,x[2]-x[1], 3137-75), nrow = 2), alternative='greater')$estimate), type='exons'))
# 
# x=read.table("SNP.count.promoter"); rownames(x)=x[,1]
# df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
# results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],46,x[2]-x[1], 3137-46), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],46,x[2]-x[1], 3137-46), nrow = 2), alternative='greater')$estimate), type='promoter'))
# 
# x=read.table("SNP.count.random"); rownames(x)=x[,1]
# df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
# results=rbind(results, cbind(Disease_or_Trait=rownames(df),df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],51,x[2]-x[1], 3137-51), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],51,x[2]-x[1], 3137-51), nrow = 2), alternative='greater')$estimate), type='random'))
# 
# results=subset(results, pvalue<0.001)
# results = results[with(results, order(type, -OR)), ]
# table(results$type)
# # HiTNEs    exons promoter   random 
# #     28       67       65       22 
# results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
# write.table(results, "eRNA.SNP.enrichments.xls", sep="\t", col.names = T, row.names = F)