#!/bin/bash
i=$1
awk '{OFS="\t"; print $1,$2,$3,$1"_"$2"_"$3}' eRNA.bed | bigWigAverageOverBed $i stdin stdout | sed 's/_/\t/g' | cut -f1-4,7,8 | intersectBed -a stdin -b gencode.v19.annotation.gtf.exons-introns.bed12.RNAseq.bigwig.bed -wao | awk '{OFS="\t"; print $0, ($17==".")?-1:(0+$17/$16); }' | groupBy -g 1,2,3,4,5,6 -c 22 -o max > $i.eRNA.measured