###########################################
## bash script for running paired-end RNAseq
###########################################
#!/bin/bash

###########################################
############## 1. Configuring
###########################################

R1=$1  # full path of R1 
R2=$2  # full path of R2 (for paired-end reads)

samplename=${R1/[.|_]R1*/}
cpu=8
index=hg19
adaptorfile=adaptor.fa
ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
Annotation_GTF=$ANNOTATION/gencode.v13.annotation.gtf
Mask_GTF=$ANNOTATION/chrM.rRNA.tRNA.gtf
BOWTIE_INDEXES=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex

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
#[ -d ../filtered ] || mkdir ../filtered
#fastq-mcf -o ../filtered/$R1 -o ../filtered/$R2 -l 16 -q 15 -w 4 -x 10 -u -P $fastqmcf $adaptorfile $R1 $R2

#cd ../filtered

#############################################
#echo "################ 3. QC"
############################################

#fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2
#rm $outputdir/$samplename/*fastqc.zip

## TOADD: RNA-seqc  (by Bin)
## 
############################################

#IFS="/" read -a tarray <<< "$R1"
#
### get the sample name
#tlen=${#tarray[@]}
#
#tsampleName=${tarray[$tlen-1]}
#
#tname=$(echo $tsampleName | awk '{gsub(/_R1.fastq/,"",$tsampleName); print $1}');
#
### RNA-SeQC output directory
#rseqDir="$outputSeQC_dir/$tname"
#mkdir $rseqDir
#
### direct the results to the directory $tmpName
#tmpName="$rseqDir/$tname"
#
#bwa aln /PHShome/bz016/neurogen/rnaSeqData/RNA-SeQCoutput/human_all_rRNA.fasta $R1 > $tmpName.sai
#bwa aln /PHShome/bz016/neurogen/rnaSeqData/RNA-SeQCoutput/human_all_rRNA.fasta $R2 > $tmpName.sai
#bwa sampe /PHShome/bz016/neurogen/rnaSeqData/RNA-SeQCoutput/human_all_rRNA.fasta $tmpName.sai $tmpName.sai $R1 $R2 > $tmpName.sam
#
#samtools view -S -b -f 4 -F 264 $tmpName.sam > $tmpName.1.bam
#samtools view -S -b -f 8 -F 260 $tmpName.sam > $tmpName.2.bam
#samtools view -S -b -f 12 -F 256 $tmpName.sam > $tmpName.3.bam
#
#java -Xmx64g -jar /PHShome/bz016/neurogen/local/picard/1.538/bin/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT I=$tmpName.1.bam I=$tmpName.2.bam I=$tmpName.3.bam O=$tmpName.strip.bam
#java -Xmx64g -jar /PHShome/bz016/neurogen/local/picard/1.538/bin/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT I=$tmpName.strip.bam F=$tmpName.strip.R1.fastq F2=$tmpName.strip.R2.fastq
#
## sequence sample name: directory + sampleName
#tmpTophatOutput="$rseqDir/TophatOutput"
#
#tophat2 --num-threads 4 --mate-inner-dist 250 --mate-std-dev 60 --library-type fr-unstranded --output-dir $tmpTophatOutput --no-novel-junc --GTF /PHShome/bz016/neurogen/local/referenceGenome/gencodeV13/gencode.v13.annotation.ExonCds.gtf /PHShome/bz016/neurogen/local/referenceGenome/hg19bt2/hg19 $tmpName.strip.R1.fastq $tmpName.strip.R2.fastq
#
#java -Xmx4g -jar /PHShome/bz016/neurogen/local/picard/1.538/bin/AddOrReplaceReadGroups.jar I=$tmpTophatOutput/accepted_hits.bam O=$tmpName.tophat.RG.bam LB=$tname PL=$tname PU=$tname SM=$tname
#
#java -jar /PHShome/bz016/neurogen/local/picard/1.538/bin/ReorderSam.jar I=$tmpName.tophat.RG.bam O=$tmpName.tophat.RG.reorder.bam R=/PHShome/bz016/neurogen/local/referenceGenome/hg19bt2/hg19.fa
#
#samtools index $tmpName.tophat.RG.reorder.bam
#
#java -jar -Xmx4g /PHShome/bz016/neurogen/local/picard/1.538/bin/MarkDuplicates.jar I=$tmpName.tophat.RG.reorder.bam O=$tmpName.tophat.RG.reorder.dup.bam METRICS_FILE=$tmpName.tophat.RG.reorder.dup.all.info.txt
#
#samtools index $tmpName.tophat.RG.reorder.dup.bam
############################################


############################################
echo "############### 4. mapping to the genome"
############################################
## tophat (output accepted_hits.sam, allow up to 100 multiple hits)
## TODO: 1) use offrated index genome_offrate3; 2)RG using HD/HC/PD etc, RG-sample use samplename
tophat -o $outputdir/$samplename --no-convert-bam --rg-id $samplename -p $cpu --read-mismatches $mm $tophat $PE_option $strand_option --max-multihits 100 --no-coverage-search genome $R1 $R2

###########################################
echo "############### 5. post-processing, format converting"
###########################################

cd $outputdir/$samplename

## sam -> bam -> sorted -> index
samtools view -Sbut $BOWTIE_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
mv accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.bam


###########################################
echo "################# 6. assembly and quantification"
###########################################

cd $outputdir/$samplename

echo "## run cufflinks to get FPKM"
cufflinks -q --no-update-check $strandoption -o ./ -p $cpu -G $Annotation_GTF -M $Mask_GTF -b $BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam

echo "## run cufflinks to do de-novo discovery"
cufflinks -q --no-update-check $strandoption -o ./denovo_cufflinks -p $cpu -g $Annotation_GTF -M $Mask_GTF -b $BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam
##echo "## run trinity to do de-novo discovery"
#Trinity.pl --output denovo_trinity --seqType fq --JM 100G --left $R1 --right $R2 --CPU $cpu
#echo "## run STAR to do de-novo discovery"
## TODO: STAR 


#echo "## run htseq for reads count"
htseq-count -m intersection-strict -t exon -i gene_id -s no -q accepted_hits.sam $ANNOTATION/gencode.v14.annotation.gtf > $samplename.hgseqcount.by.gene.tab 2> $samplename.hgseqcount.by.gene.tab.stderr
#echo "## run bedtools for reads count"
#bedtools multicov -D -split -bams accepted_hits.bam -bed $ANNOTATION/gencode.v14.annotation.bed15 > $samplename.bedtools.by.trans.tab

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
bamToBed -i $samplename.accepted_hits.bam -split > $samplename.accepted_hits.bed
sort -k1,1 $samplename.accepted_hits.bed | bedItemOverlapCount $index -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $samplename.accepted_hits.bedGraph
bedGraphToBigWig $samplename.accepted_hits.bedGraph $ANNOTATION/ChromInfo.txt $samplename.accepted_hits.bw
rm $samplename.accepted_hits.bed $samplename.accepted_hits.bedGraph

ln -s $outputdir/$samplename/isoforms.fpkm_tracking $samplename.isoforms.fpkm_tracking
ln -s $outputdir/$samplename/genes.fpkm_tracking $samplename.genes.fpkm_tracking

# gtf of assembly
echo "track name=$samplename description=$samplename visibility=pack colorByStrand='200,100,0 0,100,200'" > $samplename.transcripts.gtf
cat $outputdir/$samplename/transcripts.gtf >> $samplename.transcripts.gtf
gzip -f $samplename.transcripts.gtf

echo "!! _RNAseq.lsf job for sample $samplename is done !!"