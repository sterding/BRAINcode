###########################################
## bash script for running cuffdiff for a list of sam and gtf file
###########################################
#!/bin/sh

gtflist=$1
samlist=$2
labels=$3

# Disease (e.g.HD) vs. Control
HD=`echo $samlist | sed 's/;/\n/g' | grep "_H_" | tr '\n' ',' | sed 's/,$//g'`
Ct=`echo $samlist | sed 's/;/\n/g' | grep "_C_" | tr '\n' ',' | sed 's/,$//g'`

export BOWTIE2_INDEXES=$GENOME/$index/Sequence/Bowtie2Index/
export ANNOTATION=$GENOME/$index/Annotation/Genes

### option1: use combined gtf
#cuffcompare
cuffcompare -s $GENOME/hg19/Sequence/Chromosomes -r $ANNOTATION/genes.gtf `echo $gtflist | sed 's/;/ /g'`
# HD vs. Ct
cuffdiff -p 8 -v -L HD,Ct $strandoption -M $ANNOTATION/chrM.rRNA.tRNA.gtf -b $BOWTIE2_INDEXES/genome.fa -u -o ./cuffdiff2 cuffcmp.combined.gtf $HD $Ct &
## all vs. all
#cuffdiff -p 8 -v -L $labels $strandoption -M $ANNOTATION/chrM.rRNA.tRNA.gtf -b $BOWTIE2_INDEXES/genome.fa -u cuffcmp.combined.gtf $samlist

### option2: use the annotated gtf
# HD vs. Ct
cuffdiff -p 8 -v -L HD,Ct $strandoption -M $ANNOTATION/chrM.rRNA.tRNA.gtf -b $BOWTIE2_INDEXES/genome.fa -u -o ./cuffdiff3 $ANNOTATION/genes.gtf $HD $Ct &
## all vs. all
#cuffdiff -p 8 -v -L $labels $strandoption -M $ANNOTATION/chrM.rRNA.tRNA.gtf -b $BOWTIE2_INDEXES/genome.fa -u $ANNOTATION/genes.gtf `echo $gtflist | sed 's/;/ /g'`
