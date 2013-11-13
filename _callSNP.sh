## script to detect SNP/indel from sam/bam output of RNAseq mapping
# Authos: Xianjun Dong
# Date: 2013-11-11
# Usage: ./_callSNP.sh test.sam

input_sam=$1

# split large file into split
size=`ls -sH $input_sam | cut -f1 -d' '`
if [ "$size" -gt 10000000 ] # if larger than 10G
then
    #split -b 500M  $input_sam
    awk -f _sam2variation.awk $input_sam > ${input_sam/sam/snp}
else
    awk -f _sam2variation.awk $input_sam > ${input_sam/sam/snp}
fi

# check the relative postion of SNP on the read
textHistogram -col=4 -maxBinCount=100 ${input_sam/sam/snp} > ${input_sam/sam/relpos.hist}

# exclude SNPs located in the [0-10] region of two ends (due to the low quality)
cut -f1-2 ${input_sam/sam/snp} | uniq -c | sed 's/^ *//;s/ /\t/g' | sort -k1,1nr > ${input_sam/sam/snp.sorted}
    
# histogram of sequence depth of SNP loci
textHistogram -col=1 ${input_sam/sam/snp.sorted}


    
