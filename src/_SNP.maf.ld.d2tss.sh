#!/bin/bash

bed=$1 # input bed file
i=$2 # chr
prefix=$3

GENEpos=$GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed

fgrep -w "chr"$i $bed | intersectBed -a <(awk '{OFS="\t"; print "chr"$1,$4-1,$4,$2}' $GENOME/Annotation/Variation/1000G/Chr$i.bim) -b - -u > $prefix$i.bed
cut -f4 $prefix$i.bed > $prefix$i
# MAF
plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --freq --extract $prefix$i --silent --out $prefix$i
# number of partners with LD r2>0.8
plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --show-tags $prefix$i --list-all --tag-kb 250 --tag-r2 0.8 --silent --out $prefix$i
# distance to nearest TSS
closestBed -a $prefix$i.bed -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $GENEpos | sortBed) -D b -t first | awk '{OFS="\t"; print $4,($11<0)?-$11:$11;}' > $prefix$i.d2TSS
# combine
join -j 1 <(awk '{OFS="\t"; if(NR>1) print $2,$5}' $prefix$i.frq | sort -k1,1) <(awk '{OFS="\t"; if(NR>1) print $1,$4}' $prefix$i.tags.list | sort -k1,1) | join -j 1 - <(sort -k1,1 $prefix$i.d2TSS) | sed 's/ /\t/g' > $prefix$i.combined