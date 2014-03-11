###########################################
# bash script for running paired-end RNAseq
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 1.0
# Note: call this script in the folder of fastq file
###########################################
#!/bin/bash

###########################################
echo "["`date`"] STEP 1. Configuring"
###########################################

modulename=`basename $0`
set +o posix  #  enables the execution of process substitution e.g. http://www.linuxjournal.com/content/shell-process-redirection

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
# adaptorfile=/data/neurogen/referenceGenome/adaptor_core.fa
ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
Annotation_GTF=$ANNOTATION/gencode.v13.annotation.gtf
exons=$ANNOTATION/gencode.v13.annotation.gtf.exons.bed
Mask_GTF=$ANNOTATION/chrM.rRNA.tRNA.gtf

export BOWTIE2_INDEXES=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/
pipeline_path=$HOME/neurogen/pipeline/RNAseq/
export PATH=$pipeline_path/modules:$pipeline_path/bin:$PATH

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

###########################################
echo "["`date`"] STEP  2. quality filter: adaptor removal/clip"
###########################################

##### adaptor removal
[ -d $inputdir/../filtered ] || mkdir $inputdir/../filtered

[ -e $inputdir/../filtered/adaptor.fa ] || echo -e ">Truseq_complementary_part\nAGATCGGAAGAGC" > $inputdir/../filtered/adaptor.fa

[ ! -f $outputdir/$samplename/.status.$modulename.adaptorremoval ] && \
fastq-mcf -o $inputdir/../filtered/$R1 -o $inputdir/../filtered/$R2 -x 10 -l 15 -w 4 -u $inputdir/../filtered/adaptor.fa <(zcat $R1) <(zcat $R2) && \
touch $outputdir/$samplename/.status.$modulename.adaptorremoval 

cd $inputdir/../filtered

#############################################
echo "["`date`"] STEP 3. QC"
############################################

[ ! -f $outputdir/$samplename/.status.$modulename.fastqc ] && \
fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2 && \
rm $outputdir/$samplename/*fastqc.zip && \
touch $outputdir/$samplename/.status.$modulename.fastqc 

############################################
echo "["`date`"] STEP 4. mapping"
############################################
## tophat (output accepted_hits.sam, allow up to 100 multiple hits)
## TODO: 1) use offrated index genome_offrate3;
## Note: GATK requires @RG group has fields of ID, SM, PL, LB, PU

[ ! -f $outputdir/$samplename/.status.$modulename.mapping ] && \
tophat -o $outputdir/$samplename --no-convert-bam --rg-id $samplename --rg-sample $samplename --rg-platform ILLUMINA --rg-library $samplename --rg-platform-unit $samplename --keep-fasta-order -p $cpu --read-mismatches $mm $tophat $PE_option $strand_option --max-multihits 100 --no-coverage-search genome $R1 $R2 && \
touch $outputdir/$samplename/.status.$modulename.mapping

#[ ! -f $outputdir/$samplename/.status.$modulename.smallRNAmapping ] && \
#tophat -o $outputdir/$samplename --no-convert-bam --rg-id $samplename --rg-sample $samplename --rg-platform ILLUMINA --rg-library $samplename --rg-platform-unit $samplename --keep-fasta-order -p $cpu --read-mismatches $mm $tophat $PE_option $strand_option --max-multihits 100 --no-coverage-search genome $R1 $R2 && \
#touch $outputdir/$samplename/.status.$modulename.smallRNAmapping

###########################################
echo "["`date`"] STEP 5. post-processing, format converting"
###########################################

cd $outputdir/$samplename

[ ! -f .status.$modulename.sam2bam ] && \
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted && \
mv accepted_hits.sorted.bam accepted_hits.bam && \
samtools index accepted_hits.bam && \
touch .status.$modulename.sam2bam

[ ! -f .status.$modulename.bam2stat ] && \
samtools flagstat accepted_hits.bam > accepted_hits.bam.stat && \
touch .status.$modulename.bam2stat

### bigwig for UCSC
echo "## generating bigwig files for UCSC display"

#bamToBed -i $samplename.accepted_hits.bam -bed12 | awk '{if($1!~/_/)print}' > $samplename.accepted_hits.bed ## Note: may take more time in converting bam to sam
[ ! -f .status.$modulename.sam2bw ] && \
sam2bed -v bed12=T -v sCol=NH accepted_hits.sam | awk '{if($1!~/_/)print}' > accepted_hits.bed && \
sort -k1,1 accepted_hits.bed | bedItemOverlapCount $index -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g'> accepted_hits.bedGraph && \
bedGraphToBigWig accepted_hits.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.bw && \
touch .status.$modulename.sam2bw 

# normalized bigwig (rpm)
[ ! -f .status.$modulename.sam2normalizedbw ] && \
total_mapped_reads=`cat logs/bowtie.*t_kept_reads.log | grep -P "1\stime" | awk '{s=$1+s;}END{print s}'` && \
echo "total_mapped_reads:$total_mapped_reads" && \
awk -v tmr=$total_mapped_reads 'BEGIN{print "# total_mapped_reads="tmr;}{$4=$4*1e6/tmr; print}' accepted_hits.bedGraph > accepted_hits.normalized.bedGraph && \
bedGraphToBigWig accepted_hits.normalized.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.normalized.bw && \
touch .status.$modulename.sam2normalizedbw

# % of genome coverage vs. RNAseq density (e.g. >1rpm)
[ ! -f .status.$modulename.rpm_vs_coverage ] && \
awk 'BEGIN{max=100; UNIT=0.01; OFS="\t";}{if($0~/^#/) {print; next;} i=int($4/UNIT);if(i>max) i=max; rpm[i]+=($3-$2);}END{for(x=max;x>=0;x--) print x*UNIT, rpm[x]?rpm[x]:0;}' accepted_hits.normalized.bedGraph > accepted_hits.normalized.rpm_vs_coverage.txt && \
touch .status.$modulename.rpm_vs_coverage

# eRNA (non-generic regions with depth >10rpm and length >50bp)
#1. intersect bedGraph with all exons
#2. filter the remained region with depth cutoff
#3. filter the remained region with length cutoff
[ ! -f .status.$modulename.eRNA ] && \
intersectBed -a accepted_hits.normalized.bedGraph -b $exons -v | awk '$4>=1' | sortBed | uniq | mergeBed | awk '($3-$2)>=50' > accepted_hits.normalized.eRNA.bed && \
touch .status.$modulename.eRNA

#rm accepted_hits.bed accepted_hits.*bedGraph

###########################################
echo "["`date`"] STEP 6. call variation"
###########################################

cd $outputdir/$samplename

[ ! -f .status.$modulename.callSNP ] && \
_callSNP.sh accepted_hits.sam && \
touch .status.$modulename.callSNP

# shuilin's GATK
#[ ! -f .status.$modulename.callSNP_GATK ] && \
#_bam2vcf.sh accepted_hits.sam && \
#touch .status.$modulename.callSNP_GATK

###########################################
echo "["`date`"] STEP 7. assembly and quantification"
###########################################

cd $outputdir/$samplename

#echo "## run cufflinks to assembly (including do de-novo discovery)"
#cufflinks --no-update-check $strandoption -o ./denovo -p $cpu -g $Annotation_GTF -M $Mask_GTF accepted_hits.bam
##echo "## run trinity to do de-novo discovery"
#Trinity.pl --output denovo --seqType fq --JM 100G --left $R1 --right $R2 --CPU $cpu
#echo "## run STAR to do de-novo discovery"
## TODO: STAR

echo "## run cufflinks to get FPKM"
# Using gtf from deno assembly
# Note: "-b" option (for bias correction) can lead to segementation fault error.
[ ! -f .status.$modulename.cufflinks ] && \
cufflinks --no-update-check $strandoption -o ./ -p $cpu -G $Annotation_GTF -M $Mask_GTF --compatible-hits-norm --multi-read-correct accepted_hits.bam && \
touch .status.$modulename.cufflinks 

#echo "## run cufflinks without -M option"
#cufflinks -q --no-update-check $strandoption -o ./cufflink_wo_M -p $cpu -G $Annotation_GTF -b $BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam

echo "## run htseq for reads count"
[ ! -f .status.$modulename.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q accepted_hits.sam $Annotation_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch .status.$modulename.htseqcount

############################################
echo "["`date`"] STEP 8. for uniq"
############################################

## extract the unique mapper only
[ -d $outputdir/$samplename/uniq ] || mkdir $outputdir/$samplename/uniq
cd $outputdir/$samplename/uniq

[ ! -f $outputdir/$samplename/.status.$modulename.uniq ] && \
samtools view -SH $outputdir/$samplename/accepted_hits.sam > accepted_hits.sam && \
fgrep -w NH:i:1 $outputdir/$samplename/accepted_hits.sam >> accepted_hits.sam && \
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted && \
mv accepted_hits.sorted.bam accepted_hits.bam && \
samtools index accepted_hits.bam && \
cufflinks --no-update-check $strandoption -o ./ -p $cpu -G $Annotation_GTF -M $Mask_GTF --compatible-hits-norm accepted_hits.bam && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q accepted_hits.sam $Annotation_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
sam2bed -v bed12=T -v sCol=NH accepted_hits.sam | awk '{if($1!~/_/)print}' > accepted_hits.bed && \
sort -k1,1 accepted_hits.bed | bedItemOverlapCount $index -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > accepted_hits.bedGraph && \
mv $samplename.accepted_hits.bedGraph accepted_hits.bedGraph && \
bedGraphToBigWig accepted_hits.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.bw && \
total_mapped_reads=`wc -l accepted_hits.bed | cut -f1 -d' '` && \
echo "total_mapped_reads:$total_mapped_reads" && \
awk -v tmr=$total_mapped_reads 'BEGIN{print "# total_mapped_reads_in_million="tmr;}{$4=$4*1e6/tmr; print}' accepted_hits.bedGraph > accepted_hits.normalized.bedGraph && \
bedGraphToBigWig accepted_hits.normalized.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.normalized.bw && \
touch $outputdir/$samplename/.status.$modulename.uniq

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.callSNP ] && \
_callSNP.sh accepted_hits.sam && \
touch $outputdir/$samplename/.status.$modulename.uniq.callSNP

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.rpm_vs_coverage ] && \
awk 'BEGIN{max=100; UNIT=0.01; OFS="\t";}{if($0~/^#/) {print; next;} i=int($4/UNIT);if(i>max) i=max; rpm[i]+=($3-$2);}END{for(x=max;x>=0;x--) print x*UNIT, rpm[x]?rpm[x]:0;}' accepted_hits.normalized.bedGraph > accepted_hits.normalized.rpm_vs_coverage.txt && \
touch $outputdir/$samplename/.status.$modulename.uniq.rpm_vs_coverage

[ ! -f .status.$modulename.uniq.eRNA ] && \
intersectBed -a accepted_hits.normalized.bedGraph -b $exons -v | awk '$4>=1' | sortBed | uniq | mergeBed | awk '($3-$2)>=50' > accepted_hits.normalized.eRNA.bed && \
touch .status.$modulename.uniq.eRNA

#rm accepted_hits.bed accepted_hits.*bedGraph

############################################
echo "["`date`"] STEP 9. prepare for tracks files to display on UCSC / IGV"
############################################
#
[ -d $inputdir/../for_display ] || mkdir $inputdir/../for_display
cd $inputdir/../for_display

### others
[ ! -f $outputdir/$samplename/.status.$modulename.makelinks ] && \
# make soft link
# multi
ln -fs $outputdir/$samplename/accepted_hits.bam $samplename.multi.accepted_hits.bam && \
ln -fs $outputdir/$samplename/accepted_hits.bam.bai $samplename.multi.accepted_hits.bam.bai && \
ln -fs $outputdir/$samplename/isoforms.fpkm_tracking $samplename.multi.isoforms.fpkm_tracking && \
ln -fs $outputdir/$samplename/genes.fpkm_tracking $samplename.multi.genes.fpkm_tracking && \
ln -fs $outputdir/$samplename/hgseqcount.by.gene.tab $samplename.multi.hgseqcount.by.gene.tab && \
ln -fs $outputdir/$samplename/accepted_hits.bw $samplename.multi.accepted_hits.bw && \
ln -fs $outputdir/$samplename/accepted_hits.normalized.bw $samplename.multi.accepted_hits.normalized.bw && \
# uniq
ln -fs $outputdir/$samplename/uniq/accepted_hits.bam $samplename.uniq.accepted_hits.bam && \
ln -fs $outputdir/$samplename/uniq/accepted_hits.bam.bai $samplename.uniq.accepted_hits.bam.bai && \
ln -fs $outputdir/$samplename/uniq/isoforms.fpkm_tracking $samplename.uniq.isoforms.fpkm_tracking && \
ln -fs $outputdir/$samplename/uniq/genes.fpkm_tracking $samplename.uniq.genes.fpkm_tracking && \
ln -fs $outputdir/$samplename/uniq/hgseqcount.by.gene.tab $samplename.uniq.hgseqcount.by.gene.tab && \
ln -fs $outputdir/$samplename/uniq/accepted_hits.bw $samplename.uniq.accepted_hits.bw && \
ln -fs $outputdir/$samplename/uniq/accepted_hits.normalized.bw $samplename.uniq.accepted_hits.normalized.bw && \
# gtf of assembly
echo "track name=$samplename.multi.gtf description=$samplename.multi.gtf visibility=pack colorByStrand='200,100,0 0,100,200'" > $samplename.multi.transcripts.gtf && \
cat $outputdir/$samplename/transcripts.gtf >> $samplename.multi.transcripts.gtf && \
gzip -f $samplename.multi.transcripts.gtf && \
# for uniq
echo "track name=$samplename.uniq.gtf description=$samplename.uniq.gtf visibility=pack colorByStrand='200,100,0 0,100,200'" > $samplename.uniq.transcripts.gtf && \
cat $outputdir/$samplename/uniq/transcripts.gtf >> $samplename.uniq.transcripts.gtf && \
gzip -f $samplename.uniq.transcripts.gtf && \
## QC
ln -fs $outputdir/$samplename/$samplename.R1_fastqc $samplename.R1_fastqc && \
ln -fs $outputdir/$samplename/$samplename.R2_fastqc $samplename.R2_fastqc && \

touch $outputdir/$samplename/.status.$modulename.makelinks

echo "["`date`"] DONE: $modulename job for sample $samplename is done !!"