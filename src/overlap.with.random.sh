# script to calculate overlap for a random file (with matched length distribution as input set)
#!/bin/bash

ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
# randomly sampling
tmp=`mktemp`
bedtools shuffle -excl ~/eRNAseq/externalData/toExclude.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp
# FANTOMF CAGE
intersectBed -a $tmp -b ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed -u | wc -l >> randomoverlap.CAGE.txt
    
tmp=`mktemp`
cat ~/eRNAseq/externalData/toExclude.bed ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-3 | sortBed | mergeBed -i - | bedtools shuffle -excl - -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp

# TFBS hotspot
awk '$6>=1' ~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -c | awk '$5>=5' | wc -l >> randomoverlap.TFBS.txt
# P300
awk '$4=="EP300"' ~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.P300.txt
# Roadmpa Histone
intersectBed -a $tmp -b ~/eRNAseq/externalData/Segment/15_coreMarks_segments.E6E7E12.bed -u | wc -l >> randomoverlap.Histone.txt
# VISTA positive
grep -v NULL ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.VISTA.txt
# HCNE
intersectBed -a $tmp -b ~/eRNAseq/externalData/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | wc -l >> randomoverlap.HCNE.txt
# DNase peak from Roadmap
intersectBed -a $tmp -b ~/eRNAseq/externalData/DNase/merged.DNase.pval.signal.peaks -u | wc -l >> randomoverlap.DNase.txt
# GWAS SNPs
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | cut -f1-3 | sort -u | intersectBed -a $tmp -b stdin -u | wc -l >> randomoverlap.GWAS.txt