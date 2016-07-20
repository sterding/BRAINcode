####################################
## script to detect SNP/indel from sam/bam output of RNAseq mapping
# Authos: Xianjun Dong
# Date: 2013-11-11
# Require: sort-bed (BEDOPT), textHistogram (from Jim Kent's utility), ParaFly (only for large file)
# Usage: ./_callSNP.sh test.sam
####################################
#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` <input.sam>"
  exit
fi

pipeline_path=$HOME/neurogen/pipeline/RNAseq/
export PATH=$pipeline_path/modules:$pipeline_path/bin:$PATH

GENOME=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19

dbSNP=$GENOME/Annotation/Variation/snp137.bed.groupped.SNP
rRNA=$GENOME/Annotation/Genes/rRNA.bed
LINE=$GENOME/Annotation/Genes/LINE.bed
exons=$GENOME/Annotation/Genes/gencode.v13.annotation.gtf.exons.bed
introns=$GENOME/Annotation/Genes/gencode.v13.annotation.gtf.introns.bed

input_sam=$1

# split large file into split
size=`ls -sH $input_sam | cut -f1 -d' '`
echo "filesize is $size";
if [ "$size" -gt 10000000 ] # if larger than 10G
then
    split -b 1000M  $input_sam tmp_sampiece_
    > .paraFile;
    for i in tmp_sampiece_*; do echo "awk -f $pipeline_path/modules/_sam2variation.awk $i > tmp_snp.$i" >> .paraFile; done
    rm -f .paraFile.completed
    ParaFly -c .paraFile -CPU 8
    # merge all snp pieces
    # Note: chimeric alignment may be included, which means one read may occur multiple times in the same position. For that case, we can only count once per SNP per read. 
    cat tmp_snp* | sort -u > ${input_sam/sam/snp}
    rm tmp_snp* tmp_sampiece_*
else
    awk -f $pipeline_path/modules/_sam2variation.awk $input_sam | sort -u > ${input_sam/sam/snp}
fi

# check the relative postion of SNP on the read
textHistogram -col=4 -maxBinCount=100 ${input_sam/sam/snp} > ${input_sam/sam/snp}.relpos.hist


# get read length
readslength=`grep -v -m 10 "^@" $input_sam | awk '{print length($10)}' | sort -nr | head -n1`

# exclude SNPs located in the [0-10) region of two ends (due to the low quality)  --- a "brute force" method; GATK using a statistic model to handle this
# NOTE: even though the sam is "sorted" by Tophat, but actualy Mt hits seems not sorted. Need to double check
awk -vRL=$readslength '$4>=10 && $4<=(RL-10)' ${input_sam/sam/snp} | cut -f1-2 | sed 's/:/\t/;s/-/\t/' | bedtools groupby -g 1-4 -c 4 -o count | sort -k5,5nr > ${input_sam/sam/snp.depth}
    
# histogram of sequence depth of SNP loci
textHistogram -col=5 -minVal=10 -maxBinCount=100 -binSize=10 ${input_sam/sam/snp.depth} > ${input_sam/sam/snp.depth.hist}

# exclude SNP with <16 depth and annotate the rest
awk '{if($5>15) print; else exit;}' ${input_sam/sam/snp}.depth | sortBed >  ${input_sam/sam/snp}.depth_gt_15

# annotate SNP
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $dbSNP   -sorted -wao | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse > ${input_sam/sam/snp}_dbSNP
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $rRNA    -sorted -wao | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse > ${input_sam/sam/snp}_rRNA
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $LINE    -sorted -wao | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse > ${input_sam/sam/snp}_LINE
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $exons   -sorted -wao | sed 's/___/\t/g' | cut -f1-5,9 | uniq | groupBy -g 1,2,3,4,5 -c 6 -o collapse > ${input_sam/sam/snp}_exon
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $introns -sorted -wao | cut -f1-5,9 | uniq | groupBy -g 1,2,3,4,5 -c 6 -o collapse > ${input_sam/sam/snp}_intron

# in order of: dbSNP, exon, intron, LINE, rRNA
paste ${input_sam/sam/snp}_* | awk '{OFS="\t"; printf "%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5; for(i=6;i<=NF;i=i+6) printf "\t%s", $i; printf "\n";}' > ${input_sam/sam/snp}.annotation


exit









## NOTE:  old version
intersectBed -a ${input_sam/sam/snp}.depth_gt_15 -b $dbSNP -sorted -wao | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse | awk '{OFS="\t"; print $0, "known"}' > ${input_sam/sam/snp_known}
cut -f1-5 ${input_sam/sam/snp_known} | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | intersectBed -a stdin -b $rRNA -wo | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse | sort -k5,5nr | awk '{OFS="\t"; print $0, "rRNA"}' > ${input_sam/sam/snp_rRNA}
cut -f1-5 ${input_sam/sam/snp_known} | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | intersectBed -a stdin -b $LINE -wo | cut -f1-5,9 | groupBy -g 1,2,3,4,5 -c 6 -o collapse | sort -k5,5nr | awk '{OFS="\t"; print $0, "LINE"}' > ${input_sam/sam/snp_LINE}
# exons
cat ${input_sam/sam/snp_known} ${input_sam/sam/snp_rRNA} ${input_sam/sam/snp_LINE} | cut -f1-5 | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | intersectBed -a stdin -b $exons -wo | cut -f1-5,11,12 | sort -u | groupBy -g 1,2,3,4,5 -c 7,6 -o collapse,collapse | sort -k5,5nr | awk '{OFS="\t"; $7="exon_"$7; print}' > ${input_sam/sam/snp_exon}
# intronic
cat ${input_sam/sam/snp_known} ${input_sam/sam/snp_rRNA} ${input_sam/sam/snp_LINE} | cut -f1-5 | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | intersectBed -a stdin -b $exons -wo | cut -f1-5,11,12 | sort -u | groupBy -g 1,2,3,4,5 -c 7,6 -o collapse,collapse | sort -k5,5nr | awk '{OFS="\t"; $7="exon_"$7; print}' > ${input_sam/sam/snp_exon}
# proximal (1M)
cat ${input_sam/sam/snp_known} ${input_sam/sam/snp_rRNA} ${input_sam/sam/snp_LINE} | cut -f1-5 | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | intersectBed -a stdin -b $exons -wo | cut -f1-5,11,12 | sort -u | groupBy -g 1,2,3,4,5 -c 7,6 -o collapse,collapse | sort -k5,5nr | awk '{OFS="\t"; $7="exon_"$7; print}' > ${input_sam/sam/snp_exon}

cat ${input_sam/sam/snp_known} ${input_sam/sam/snp_rRNA} ${input_sam/sam/snp_LINE} ${input_sam/sam/snp_exon} | cut -f1-5 | sort - ${input_sam/sam/snp}.depth_gt_15 | uniq -u | sort -k5,5nr | awk '{OFS="\t"; print $0, "others", "others"}' > ${input_sam/sam/snp_others}

cat  ${input_sam/sam/snp}_* > ${input_sam/sam/snp}.annotation
rm ${input_sam/sam/snp}_*
