# script to calculate overlap for a random file
#!/bin/bash

ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
sortBed -i $ANNOTATION/exons.meta.bed | mergeBed -i - | cat - $ANNOTATION/hg19.gap.bed  > /tmp/toexclude.bed
# randomly sampling
tmp=`mktemp`
bedtools shuffle -excl /tmp/toexclude.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp
# TFBS hotspot
awk '$6>=1' ../TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -c | awk '$5>5' | wc -l >> randomoverlap.TFBS.txt
# P300
awk '$4=="EP300"' ../TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.P300.txt
# FANTOMF CAGE
intersectBed -a $tmp -b ../CAGE/permissive_enhancers.bed -u | wc -l >> randomoverlap.CAGE.txt
# Roadmpa Histone
intersectBed -a $tmp -b ../Segment/15_coreMarks_segments.E6E7E12.bed -u | wc -l >> randomoverlap.Histone.txt
# VISTA positive
grep -v NULL /PHShome/xd010/projects/PD/results/eRNA/externalData/VISTA/hg19.tested_regions.bed | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.VISTA.txt
# HCNE
intersectBed -a $tmp -b ../Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | wc -l >> randomoverlap.HCNE.txt
# DNase enhancer from Roadmap
intersectBed -a $tmp -b ../DNase/regions_enh_merged.brain.bed -u | wc -l >> randomoverlap.DNase.txt
# GWAS SNPs
cut -f1-3 $GENOME/Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.bed | sort -u | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.GWAS.txt