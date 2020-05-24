# script to parse bam/sam file to see the porpotion of reads mapped to different annotation regions

if [ $# -ne 1 ]
then
  echo "Usage: $0 accepted_hits.bam"
  exit
fi

inputbam=$1

# annotation files (using export because of parallel 'xxx')
export ANNOTATION=$GENOME/Annotation/Genes
export inputbam=$1 

# total reads
echo "total:" `samtools view $inputbam -c`
echo "rRNA:" `samtools view $inputbam -c -L $ANNOTATION/rRNA.bed`
echo "mtRNA:" `samtools view $inputbam -c chrM`

#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.intergenic.bed
#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools flank -g $ANNOTATION/ChromInfo.txt -b 5000 | sortBed | mergeBed > $ANNOTATION/gencode.v19.annotation.intergenic.flank5k.bed
#cut -f1-3 $ANNOTATION/gencode.v19.annotation.bed12 | sortBed | mergeBed | bedtools slop -g $ANNOTATION/ChromInfo.txt -b 5000 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.intergenic.notneargenes.bed
#make soft link to intergenic.bed intergenic_neargene.bed intergenic_notneargene.bed

#echo -e "chrM\t0\t16571" | cat - $ANNOTATION/rRNA.bed | cut -f1-3 | sortBed | mergeBed | bedtools complement -g $ANNOTATION/ChromInfo.txt > $ANNOTATION/gencode.v19.annotation.non_rRNA_mt.bed
# ln -fs gencode.v19.annotation.non_rRNA_mt.bed non_rRNA_mt.bed
samtools view $inputbam -b -1 -L $ANNOTATION/non_rRNA_mt.bed -o $inputbam.non-rRNA-mt.bam
samtools index $inputbam.non-rRNA-mt.bam
echo "total_non_rRNA_mt:" `samtools view $inputbam.non-rRNA-mt.bam -c`
echo -e "intergenic\nintergenic_notneargene\nintrons\nexons\nutr5\nutr3\nLINE\nSINE" | parallel 'echo {}: `samtools view $inputbam.non-rRNA-mt.bam -c -L $ANNOTATION/{}.bed`'