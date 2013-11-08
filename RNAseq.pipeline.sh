####################################
# Pipeline for RNAseq data analysis
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 1.0
####################################
#!/bin/bash

############
## 0. setting 
############
n_CPU=8
reference_version=hg19
ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
Annotation_GTF=$ANNOTATION/gencode.v13.annotation.gtf
Mask_GTF=$ANNOTATION/chrM.rRNA.tRNA.gtf
BOWTIE_INDEXES=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex

pipeline_path=$HOME/neurogen/pipeline/RNAseq/
export PATH=$pipeline_path:$PATH

##TODO: test if the executable program are installed 
# bowtie, tophat, cufflinks, htseq-count, bedtools, samtools, RNA-seQC ... 

input_dir=$1  # $HOME/neurogen/xdong/rnaseq_PD/rawfiles
output_dir=$input_dir/../run_output
[ -d $output_dir ] || mkdir $output_dir


## Add a folder for output from running RNA-SeQC pipeline from Nathlie Broad
outputSeQC_dir=$input_dir/run_RNA-SeQC
[ -d $outputSeQC_dir ] || mkdir $outputSeQC_dir

adaptor_file=adaptor.fa

## parameters of RNAseq library
#phred
bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger"; fastqmcf="33"; trimmomatic="-phred33"
#Pair-end option
PE="--mate-inner-dist 50 --mate-std-dev 20"
#strand
strandoption="--library-type fr-unstranded"

## mapping
#mismatch
mm=2

## hpcc cluster setting
email="-u sterding.hpcc@gmail.com -N"
cpu=8
memory=40000 # unit in Kb, e.g. 20000=20G


############
## 1. QC/mapping/assembly/quantification for all samples in the input dir  (Tophat/Cufflink/Htseq-count)
############
cd $input_dir

c=0;h=0;gtflist="";samlist=""; labels=""

for i in *R1.fastq;  # for test purpose only use samples starting with 2*
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -o $output_dir/$samplename/_RNAseq.log -q big-multi -n $cpu -M $memory -R rusage[mem=$memory] $email _RNAseq.sh $R1 $R2
    
    #jobid=`bsub RNAseq.lsf $R1 $R2 | cut -f3 -d' '`
    #echo "Your job is submitted (jobID: $jobid) with SGE script at $output_dir/$samplename/$samplename.sge"

    gtflist="$gtflist;$output_dir/$samplename/transcripts.gtf"
    samlist="$samlist;$output_dir/$samplename/accepted_hits.sam"
    if [ "$labels" == "" ];
        then
            labels="$samplename";
        else
            labels="$labels,$samplename"
    fi
done

# set break point here to wait until all samples are completedly processed.
exit

############
## 2. TODO: Outlier analysis (clustering, visualization) -- by Bin
############
bsub _outlier.sh 

############
## 3. factor analysis to identify the hidden covariates (PEER)
############
bsub Rscript _factor_analysis.R

############
## 4. identify differentially expressed genes (cuffdiff and DEseq), incoperating the hidden covariates from PEER
############
[ -d $output_dir/DE_cuffdiff ] || mkdir $output_dir/DE_cuffdiff
cd $output_dir/DE_cuffdiff

bsub -o $output_dir/$samplename/_DE_cuffdiff.log -q long -n $cpu -R rusage[mem=$memory] -u $email -N _DE_cuffdiff.sh $gtflist $samlist $labels
bsub Rscript _DE_DEseq.R $output_dir PD Ct $ANNOTATION


############
## 5. eQTL (PEER)
############
## pre-requirisition: call SNP/variation ahead  -- by Shuilin
bsub _bam2vcf.sh $bamfile # by Shuilin

# eQTL
bsub Rscript _eQTL.R

############
## 6. pathway analysis (SPIA)
############
bsub Rscript _pathway_analysis.R
