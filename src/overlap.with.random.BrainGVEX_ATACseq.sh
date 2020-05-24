# script to calculate overlap for a random file (with matched length distribution as input set)
#!/bin/bash

ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
# randomly sampling
tmp=`mktemp`
bedtools shuffle -excl ~/eRNAseq/externalData/toExclude.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp
# sknsh_ENCODE
awk '$4>=5' ~/neurogen/psychENCODE/BrainGVEX/ATAC-seq/merged.Peaks.bedGraph | intersectBed -a $tmp -b - -u | wc -l >> randomoverlap.BrainGVEX_ATACseq.txt
    
