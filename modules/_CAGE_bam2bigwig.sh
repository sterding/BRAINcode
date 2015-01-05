# convert CAGE bam to bigwig
#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Script to generate bigwig files for CAGE bam input"
  echo "=================================================="
  echo "Usage: $0 input.bam"
  echo "Output: input.plus.bw and input.minus.bw"
  exit
fi

inputbam=$1

# use the 5' 1st nt as representative

#samtools view $inputbam | sam2bed -v bed12=F -vXSstrand=F | awk '{OFS="\t"; s=($6=="+")?$2:($3-1); e=($6=="+")?($2+1):$3; print $1,s,e,$4,$5,$6;}' > ${inputbam/bam/bed}

# use the 5' 2nd nt as representative

samtools view $inputbam | sam2bed -v bed12=F -vXSstrand=F | awk '{OFS="\t"; s=($6=="+")?($2+1):$3; e=($6=="+")?($2+2):($3+1); print $1,s,e,$4,$5,$6;}' > ${inputbam/bam/bed}
bam2bigwig.sh ${inputbam/bam/bed} -split