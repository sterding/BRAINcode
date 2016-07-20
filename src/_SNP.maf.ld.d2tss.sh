#!/bin/bash
# usage: bsub -J "maf[1-22]" -q normal -n 1 bash $pipeline_path/src/_SNP.maf.ld.d2tss.sh

bed=$1 # input bed file
prefix=$2
i=$LSB_JOBINDEX # chr

# debug
# bed=/data/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed; prefix='ALL_chr'; i=1

GENEpos=$GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed

fgrep -w "chr"$i $bed | intersectBed -a <(awk '{OFS="\t"; print "chr"$1,$4-1,$4,$2}' $GENOME/Annotation/Variation/1000G/Chr$i.bim) -b - -u > $prefix$i.bed
cut -f4 $prefix$i.bed > $prefix$i
# MAF
/source/plink/1.07/plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --freq --extract $prefix$i --silent --out $prefix$i
# number of partners with LD r2>0.8
/source/plink/1.07/plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --show-tags $prefix$i --list-all --tag-kb 250 --tag-r2 0.8 --silent --out $prefix$i
# distance to nearest TSS
closestBed -a $prefix$i.bed -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $GENEpos | sortBed) -D b -t first | awk '{OFS="\t"; print $4,($11<0)?-$11:$11;}' > $prefix$i.d2TSS
# combine
join -j 1 <(awk '{OFS="\t"; if(NR>1) print $2,$5}' $prefix$i.frq | sort -k1,1) <(awk '{OFS="\t"; if(NR>1) print $1,$4}' $prefix$i.tags.list | sort -k1,1) | join -j 1 - <(sort -k1,1 $prefix$i.d2TSS) | sed 's/ /\t/g' > $prefix$i.combined
# clean up
rm $prefix$i $prefix$i.bed $prefix$i.frq $prefix$i.tags.list $prefix$i.d2TSS 