### Pipeline for RNAseq data analysis

####################################
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 1.0
####################################

############
## 0. setting 
############
n_CPU=8
reference_version=hg19
ANNOTATION=/data/neurogen/commonData/gencode.v14.annotation.bed15
BOWTIE_INDEXES=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex
BIN=$HOME/projects/PD/src

input_dir=$1  # $HOME/neurogen/xdong/rnaseq_PD/rawfiles
output_dir=$input_dir/../run_output
[ -d $output_dir ] || mkdir $output_dir

adaptor_file=adaptor.fa

## parameters of RNAseq library
#phred
bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger"; fastqmcf="33"; trimmomatic="-phred33"
#Pair-end option
PE="--mate-inner-dist 250 --mate-std-dev 60"
#strand
strandoption="--library-type fr-unstranded"

## mapping
#mismatch
mm=2

############
## 1. QC/mapping/assembly/quantification for all samples in the input dir  (Tophat/Cufflink/Htseq-count)
############
cd $input_dir

c=0;h=0;gtflist="";samlist=""; labels=""

for i in *.R1.fastq;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub RNAseq.lsf $R1 $R2
    
    jobid=`qsub $output_dir/$samplename/$samplename.sge | cut -f3 -d' '`

    echo "Your job is submitted (jobID: $jobid) with SGE script at $output_dir/$samplename/$samplename.sge"

    gtflist="$gtflist;$output_dir/$samplename/transcripts.gtf"
    samlist="$samlist;$output_dir/$samplename/accepted_hits.sam"
    if [ "$labels" == "" ];
        then
            labels="$samplename";
        else
            labels="$labels,$samplename"
    fi
done


############
## 2. TODO: Outlier analysis (clustering, visualization) -- by Bin
############
bsub outlier.sh 

############
## 3. factor analysis to identify the hidden covariates (PEER)
############
bsub Rscript _factor_analysis.R

############
## 4. identify differentially expressed genes (cuffdiff and DEseq), incoperating the hidden covariates from PEER
############
[ -d $output_dir/DE_cuffdiff ] || mkdir $output_dir/DE_cuffdiff
cd $output_dir/DE_cuffdiff

bsub _DE_cuffdiff.lsf $gtflist $samlist $labels
bsub Rscript _DE_DEseq.R $output_dir PD Ct $ANNOTATION


############
## 5. eQTL (PEER)
############
## pre-requirisition: call SNP/variation ahead  -- by Shuilin
bsub Rscript _eQTL.R

############
## 6. pathway analysis (SPIA)
############
bsub Rscript _pathway_analysis.R
