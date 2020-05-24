# bash script to fix bug in "bedtools coverage -split -f -counts" (see https://github.com/arq5x/bedtools2/issues/673)

module load bedtools2/2.27.1; 

inputbam=$1
inputbed=$2
genomeindex=$3
chr=$4

TMPDIR=/data/neurogen/tmp/xd010

sampleName=`echo $inputbam | sed 's/.*output\/\(.*\)\/accepted_hits.*/\1/g'`

[ ! -f $TMPDIR/$sampleName.$chr.sorted.bed ] && \
samtools view -b -o $TMPDIR/$sampleName.$chr.bam $inputbam $chr && \
bedtools bamtobed -split -i $TMPDIR/$sampleName.$chr.bam > $TMPDIR/$sampleName.$chr.bed && \
LC_ALL=C sort -k2,2n --buffer-size=1G $TMPDIR/$sampleName.$chr.bed > $TMPDIR/$sampleName.$chr.sorted.bed

bedtools coverage -a $inputbed.$chr -b $TMPDIR/$sampleName.$chr.sorted.bed -counts -sorted -g $genomeindex -f 1.0 > $TMPDIR/$inputbed.$sampleName.$chr
