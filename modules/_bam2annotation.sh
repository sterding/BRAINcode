# script to parse bam/sam file to see the porpotion of reads mapped to different annotation regions

inputbam=$1

# annotation files
export ANNOTATION=$GENOME/Annotation/Genes
export inputbam=$1

# total reads
echo "total:" `samtools view $inputbam -c`
echo "rRNA:" `grep -v chrM $ANNOTATION/rRNA.bed | samtools view $inputbam -c -L -`
echo "mtRNA:" `samtools view $inputbam -c chrM`

#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.intergenic.bed
#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools flank -g $ANNOTATION/ChromInfo.txt -b 5000 | sortBed | mergeBed > $ANNOTATION/gencode.v19.annotation.intergenic.flank5k.bed
#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools slop -g $ANNOTATION/ChromInfo.txt -b 5000 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.intergenic.notneargenes.bed
#make soft link to intergenic.bed intergenic_neargene.bed intergenic_notneargene.bed

#echo -e "chrM\t0\t16571" | cat - $ANNOTATION/rRNA.bed | cut -f1-3 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.non_rRNA_mt.bed
samtools view $inputbam -b -1 -L $ANNOTATION/gencode.v19.annotation.non_rRNA_mt.bed -o $inputbam.non-rRNA-mt.bam
echo "total_non-rRNA-mt:" `samtools view $inputbam.non-rRNA-mt.bam -c`
echo -e "intergenic\nintergenic_notneargene\nintrons\nexons\nutr5\nutr3\nLINE\nSINE" | parallel 'echo {}: `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/{}.bed`'


#echo "intergenic:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/intergenic.bed `
#echo "intron:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/introns.bed`
#echo "exon:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/exons.bed`
#echo "utr5:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/5utr.bed`
#echo "utr3:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/3utr.bed`
#echo "LINE:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/LINE.bed`
#echo "SINE:" `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/SINE.bed`
