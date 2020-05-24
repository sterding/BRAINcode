# script to calculate overlap for a random file (with matched length distribution as input set)
#!/bin/bash

ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
# randomly sampling
tmp=`mktemp`
bedtools shuffle -excl ~/eRNAseq/externalData/toExclude.bed -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt > $tmp
# sknsh_ENCODE
zcat ~/eRNAseq/externalData/sknsh_ENCODE/sknsh.*gz | cut -f1-3 | sortBed | mergeBed | intersectBed -a $tmp -b - -u | wc -l >> randomoverlap.sknsh_ENCODE.txt
    
