#!/usr/bin/env bash
# bash script to scan a fa file with a meme motif file and parse result to indivudal MOTIF.bed file
# Usage:
# cd ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/; 
# for i in MA*.meme; do bsub -q short -n 1 scan_motif.sh $i /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa; done
# cat MA*.bed | LC_ALL=C sort --parallel=24 --buffer-size=5G -k1,1 -k2,2n > ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.in_genome.fa.bed

motif_meme=$1
seq_fasta=$2
output_file=$3

[ -z "$output_file" ] && output_file=`basename $motif_meme`"_in_"`basename $seq_fasta`".bed"

#echo $0 $motif_meme $seq_fasta $output_file

# for now we just record the chr, start, and end. More info can be included

#fimo --text $motif_meme $seq_fasta | awk '$1 ~ /MA/ {OFS="\t"; print $2,$3,$4 >> $1".bed"}'
tmp_file=`mktemp -p ~/neurogen/tmp/xd010`
echo $tmp_file
fimo --text $motif_meme $seq_fasta > $tmp_file
awk '$1 ~ /MA/ {OFS="\t"; print $2,$3,$4,$1}' $tmp_file > $output_file