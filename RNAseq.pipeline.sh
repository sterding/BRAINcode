####################################
# Pipeline for RNAseq data Analysis
# Authors: Bioinformatics Team @ Scherzer's lab 
# Email: xdong@rics.bwh.harvard.edu
# Date: 9/16/2013
# Version: 1.0
####################################
#!/bin/bash

modulename=`basename $0`
set +o posix  #  enables the execution of process substitution e.g. http://www.linuxjournal.com/content/shell-process-redirection
STEP=0

if [ $# -ne 1 ]
then
  echo "Usage: $HOME/neurogen/pipeline/RNAseq/RNAseq.pipeline.sh /data/neurogen/rnaseq_PD/rawfiles"
  exit
fi

############
## 0. setting 
############
reference_version=hg19
ANNOTATION=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes
Annotation_GTF=$ANNOTATION/gencode.v13.annotation.gtf
Mask_GTF=$ANNOTATION/chrM.rRNA.tRNA.gtf
BOWTIE_INDEXES=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex
pipeline_path=$HOME/neurogen/pipeline/RNAseq/
export PATH=$pipeline_path/modules:$pipeline_path/bin:$PATH

## hpcc cluster setting
email="-u sterding.hpcc@gmail.com -N"
cpu="-n 8"
memory="-M 4000 -R rusage[mem=4000]" # unit in Kb, e.g. 20000=20G

##TODO: test if the executable program are installed 
# bowtie, tophat, cufflinks, htseq-count, bedtools, samtools, RNA-seQC ... 

# project folders
input_dir=$1  # $HOME/neurogen/xdong/rnaseq_PD/rawfiles

# create the subfolders (e.g. filtered, run_output, for_display, results)
filtered_dir=$input_dir/../filtered
[ -d $filtered_dir ] || mkdir $filtered_dir

output_dir=$input_dir/../run_output
[ -d $output_dir ] || mkdir $output_dir

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results 
[ -d $result_dir ] || mkdir $result_dir


############
## 1. QC/mapping/assembly/quantification for all samples in the input dir  (Tophat/Cufflink/Htseq-count)
############
cd $input_dir

c=0;h=0;gtflist="";samlist=""; labels=""

for i in *_4.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q big-multi $cpu $memory $email _RNAseq.sh $R1 $R2;

    gtflist="$gtflist;$output_dir/$samplename/transcripts.gtf"
    samlist="$samlist;$output_dir/$samplename/accepted_hits.sam"
    if [ "$labels" == "" ];
        then
            labels="$samplename";
        else
            labels="$labels,$samplename"
    fi
done

exit

############
## 2. draw aggregation plot for the RNAseq density in the genetic region
# Note: see note in bigWigAverageOverBed_81bins.sh for requisition 
############
[ -d $result_dir/aggregationPlot ] || mkdir $result_dir/aggregationPlot
cd $result_dir/aggregationPlot
LOG=$result_dir/aggregationPlot/_aggPlot.log
for i in $fordisplay_dir/*uniq*ed.bw;
do
    bsub -J aggPlot -oo $LOG.$i -eo $LOG.$i -q big-multi $email -M 4000 -R rusage[mem=4000] -n 2 aggregationPlot_getBinSignal.sh $GENOME/Annotation/Genes/gencode.v19.annotation.bed12.mRNA.81bins.bed 81 $i;
done

# set break point here to wait until all samples are completedly processed.
exit

############
## 2. Added cluster procedure -- by Bin
############
bsub Rscript _clustComRNASeq.R $output_dir $result_dir

############
# RNA-SeQC
#
# make samplelist_file.txt

./makeSamplelistFile.sh "$output_dir/*/acce*.bam" > $resultOutput_dir/RNA-SEQCfolder/samplelist_file.txt

# gencodePlusMask.v13.annotation.gtf is a combination of 
# gencode.v13.annotation.karotyped.gtf and chrM.rRNA.tRNA.gtf

bsub -J runRnaSeqc -q big-multi -n 4 -R 'rusage[mem=10000]' "java -jar -Xmx64g RNA-SeQC_v1.1.7.jar -s $resultOutput_dir/RNA-SEQCfolder/samplelist_file.txt -t /PHShome/bz016/neurogen/rnaseq_PD/results/RNA-SEQCfolder/gencodePlusMask.v13.annotation.gtf -r /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa -o $resultOutput_dir/RNA-SEQCfolder/RNASeQCoutput"
############


############
## 3. factor analysis to identify the hidden covariates (PEER)
############
bsub Rscript _factor_analysis.R

############
## 4. identify differentially expressed genes (cuffdiff and DEseq), incoperating the hidden covariates from PEER
############
[ -d $result_dir/DE_cuffdiff ] || mkdir $result_dir/DE_cuffdiff
cd $result_dir/DE_cuffdiff

bsub -o _DE_cuffdiff.log -q long $cpu $memory $email _DE_cuffdiff.sh $gtflist $samlist $labels

[ -d $result_dir/DE_DESeq2 ] || mkdir $result_dir/DE_DESeq2
cd $result_dir/DE_DESeq2
# Save the covariance table from Google Doc
# wget --no-check-certificate -qO - "https://docs.google.com/spreadsheet/ccc?key=0Aumm3V3g3dF7dEFnZ2pPQjlheXlZand6YWUxeF9PMUE&gid=5&output=txt" > covariances.tab
bsub -o _DE_DESeq2.log -q long $cpu $memory $email Rscript _DE_DESeq2.R $output_dir PD Ct $ANNOTATION

############
## 5. eQTL 
############
[ -d $result_dir/eQTL ] || mkdir $result_dir/eQTL
cd $result_dir/eQTL
## pre-requirisition: call SNP/variation ahead  -- by Shuilin
bsub _bam2vcf.sh $bamfile # by Shuilin

# eQTL
bsub Rscript _eQTL.R

############
## 6. pathway analysis (SPIA)
############
bsub Rscript _pathway_analysis.R
