###########################################
## bash script for running paired-end RNAseq
## Note: call this script in the folder of fastq file
###########################################
#!/bin/bash

###########################################
############## 1. Configuring
###########################################

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` HC_BN10-39_2.R1.fastq.gz HC_BN10-39_2.R2.fastq.gz"
  exit
fi

R1=$1  # filename of R1 
R2=$2  # filename of R2 (for paired-end reads)

samplename=${R1/[.|_]R1*/}
cpu=8
index=hg19
adaptorfile=/data/neurogen/referenceGenome/adaptor.fa
ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
Annotation_GTF=$ANNOTATION/gencode.v13.annotation.gtf
Mask_GTF=$ANNOTATION/chrM.rRNA.tRNA.gtf
export BOWTIE2_INDEXES=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/

#============= mapping options
#phred
bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger"; fastqmcf="33"; trimmomatic="-phred33"
#mismatch
mm=2
#PE option
PE_option="--mate-inner-dist 50 --mate-std-dev 20"  ## Shuijin found that 50/20 can get higher mappability
#strand
strand_option="--library-type fr-unstranded"

inputdir=$PWD
outputdir=$inputdir/../run_output
[ -d $outputdiri/$samplename ] || mkdir -p $outputdir/$samplename

#ln -sf /tmp/$LSB_JOBID.out $outputdir/$samplename/_RNAseq.log

###########################################
echo "###############  2. quality filter: adaptor removal/clip"
###########################################

##### adaptor removal
[ -d $inputdir/../filtered ] || mkdir $inputdir/../filtered
fastq-mcf -o $inputdir/../filtered/$R1 -o $inputdir/../filtered/$R2 -w 4 -x 10 -u -P $fastqmcf $adaptorfile $R1 $R2

cd ../filtered

#############################################
echo "################ 3. QC"
############################################

fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2
rm $outputdir/$samplename/*fastqc.zip

############################################
echo "############### 4. mapping to the genome"
############################################
## tophat (output accepted_hits.sam, allow up to 100 multiple hits)
## TODO: 1) use offrated index genome_offrate3; 
tophat -o $outputdir/$samplename --no-convert-bam --rg-id $samplename --rg-sample $samplename --keep-fasta-order -p $cpu --read-mismatches $mm $tophat $PE_option $strand_option --max-multihits 100 --no-coverage-search genome $R1 $R2

###########################################
echo "############### 5. post-processing, format converting"
###########################################

cd $outputdir/$samplename

## sam -> bam -> sorted -> index
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
mv accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.bam


###########################################
echo "################# 6. assembly and quantification"
###########################################

cd $outputdir/$samplename

echo "## run cufflinks to assembly (including do de-novo discovery)"
cufflinks -q --no-update-check $strandoption -o ./denovo -p $cpu -g $Annotation_GTF -M $Mask_GTF accepted_hits.bam
##echo "## run trinity to do de-novo discovery"
#Trinity.pl --output denovo --seqType fq --JM 100G --left $R1 --right $R2 --CPU $cpu
#echo "## run STAR to do de-novo discovery"
## TODO: STAR

echo "## run cufflinks to get FPKM"
# Using gtf from deno assembly
# Note: "-b" option (for bias correction) can lead to segementation fault error. 
cufflinks -q --no-update-check $strandoption -o ./ -p $cpu -G denovo/transcripts.gtf -M $Mask_GTF --compatible-hits-norm --multi-read-correct accepted_hits.bam

#echo "## run cufflinks without -M option"
#cufflinks -q --no-update-check $strandoption -o ./cufflink_wo_M -p $cpu -G $Annotation_GTF -b $BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam


#echo "## run htseq for reads count"
htseq-count -m intersection-strict -t exon -i gene_id -s no -q accepted_hits.sam denovo/transcripts.gtf > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr
#echo "## run bedtools for reads count"
#bedtools multicov -D -split -bams accepted_hits.bam -bed $ANNOTATION/gencode.v14.annotation.bed15 > bedtools.by.trans.tab

############################################
echo "############### 7. prepare for tracks files to display on UCSC / IGV"
############################################
#
[ -d $inputdir/../for_display ] || mkdir $inputdir/../for_display
cd $inputdir/../for_display

## make index for the (sorted) BAM
ln -s $outputdir/$samplename/accepted_hits.bam $samplename.accepted_hits.bam
ln -s $outputdir/$samplename/accepted_hits.bam.bai $samplename.accepted_hits.bam.bai

## QC
#mv $outputdir/$samplename/*_fastqc ./

# bigwig for UCSC
echo "## generating bigwig files for UCSC display"
bamToBed -i $outputdir/$samplename/accepted_hits.sam -split > $samplename.accepted_hits.bed
sort -k1,1 $samplename.accepted_hits.bed | bedItemOverlapCount $index -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $samplename.accepted_hits.bedGraph
bedGraphToBigWig $samplename.accepted_hits.bedGraph $ANNOTATION/ChromInfo.txt $samplename.accepted_hits.bw
rm $samplename.accepted_hits.bed $samplename.accepted_hits.bedGraph

ln -s $outputdir/$samplename/isoforms.fpkm_tracking $samplename.isoforms.fpkm_tracking
ln -s $outputdir/$samplename/genes.fpkm_tracking $samplename.genes.fpkm_tracking

# gtf of assembly
echo "track name=$samplename description=$samplename visibility=pack colorByStrand='200,100,0 0,100,200'" > $samplename.transcripts.gtf
cat $outputdir/$samplename/denovo/transcripts.gtf >> $samplename.transcripts.gtf
gzip -f $samplename.transcripts.gtf

echo "!! _RNAseq.lsf job for sample $samplename is done !!"
