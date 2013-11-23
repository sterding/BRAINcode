## script to detect SNP/indel from sam/bam output of RNAseq mapping
# Authos: Xianjun Dong
# Date: 2013-11-11
# Require: sort-bed (BEDOPT), textHistogram (from Jim Kent's utility), ParaFly (only for large file)
# Usage: ./_callSNP.sh test.sam

input_sam=$1

# split large file into split
size=`ls -sH $input_sam | cut -f1 -d' '`
echo "filesize is $size";
if [ "$size" -gt 10000000 ] # if larger than 10G
then
    split -b 1000M  $input_sam tmp_sampiece_
    for i in tmp_sampiece_*; do echo "awk -f _sam2variation.awk $i > tmp_snp.$i" >> .paraFile; done
    ParaFly -c .paraFile -CPU 8
    # merge all snp pieces
    cat tmp_snp* > ${input_sam/sam/snp}
else
    awk -f _sam2variation.awk $input_sam > ${input_sam/sam/snp}
fi

# check the relative postion of SNP on the read
textHistogram -col=4 -maxBinCount=100 ${input_sam/sam/snp} > ${input_sam/sam/relpos.hist}

# exclude SNPs located in the [0-10] region of two ends (due to the low quality)
# NOTE: even though the sam is "sorted" by Tophat, but actualy Mt hits seems not sorted. Need to double check
awk '{if($4>10 && $4<40) print}' ${input_sam/sam/snp} | cut -f1-2 | sed 's/:/\t/;s/-/\t/' | sort --buffer-size=2G -k1,1 -k2,2n -k4,4 | bedtools groupby -g 1-4 -c 4 -o count | sort -k5,5nr > ${input_sam/sam/snp.depth}
    
# histogram of sequence depth of SNP loci
textHistogram -col=1 -minVal=10 -maxBinCount=1000 -binSize=10 ${input_sam/sam/snp.depth} > ${input_sam/sam/snp.depth.hist}

# exclude SNP with <10 depth and annotate the rest
rRNA=~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/rRNA.bed
LINE=~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/LINE.bed
dbSNP=~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/snp137.bed

## UNFINISHED

intersectBed -a ${input_sam/sam/snp.depth} -b $dbSNP -wo |
intersectBed -a ${input_sam/sam/snp.depth} -b $dbSNP -v | awk '{if($5>10) print; else exit;}' | intersectBed -a stdin -b $rRNA |
intersectBed -a ${input_sam/sam/snp.depth} -b $dbSNP -v | awk '{if($5>10) print; else exit;}' | intersectBed -a stdin -b $LINE | 