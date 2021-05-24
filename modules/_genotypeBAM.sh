# Script to get "genotype" for a list of SNPs using RNAseq bam alignment
# author: Dong Xianjun
# Date: 2016-05-19
# Version: 1.0
# Usage: bsub -q big -n 4 bash _genotypeBAM.sh ~/neurogen/rnaseq_PD/results/RNAvsDNA/Overlap.Exon.SNP.bed ~/neurogen/rnaseq_PD/results/RNAvsDNA/HC_SNDA.bamlist
# Input: a list of position in (chr pos) format or BED format, and a list of bam files
# Note: the intermidiate file starting with gatk_xxxx can be deleted later.

position_file=$1
bam_list_file=$2

# call vcf from bam
samtools mpileup -uBIv -l $position_file -f $GENOME/Sequence/WholeGenomeFasta/genome.fa -o $position_file.vcf -t AD,ADR,ADF,DP,INFO/AD,INFO/ADF,INFO/ADR -b $bam_list_file

# get genotype per sample (multiallelic-caller; no indel; add GQ; minimal DP per sample:15;)
bcftools call -m -V indels -f GQ -g 15 $position_file.vcf > $position_file.call.vcf

# convert to tped using plink
plink --vcf $position_file.call.vcf --recode transpose --const-fid --out $position_file

# convert to customerized tped format 
awk '{print $2}' $position_file.tfam | rowsToCols stdin stdout | awk 'BEGIN{OFS="\t"}{print "Chr","START","END","ID", $0}' > $position_file.tped.txt
awk '{printf "%s\t%s\t%s\t%s", "chr"$1, $4-1, $4, $2; for(i=5;i<=NF;i=(i+2)) { j=i+1; printf "\t%s",$i$j;} printf "\n"; }' $position_file.tped >> $position_file.tped.txt

echo "DONE!"