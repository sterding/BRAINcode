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

########################
## 0. setting 
########################
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
memory="-M 40000 -R rusage[mem=40000]" # unit in Kb, e.g. 20000=20G

##TODO: test if the executable program are installed 
# bowtie, tophat, cufflinks, htseq-count, bedtools, samtools, RNA-seQC ... 

# project folders
input_dir=$1  # input_dir=/data/neurogen/rnaseq_PD/rawfiles

# create the subfolders (e.g. filtered, run_output, for_display, results)
filtered_dir=$input_dir/../filtered
[ -d $filtered_dir ] || mkdir $filtered_dir

output_dir=$input_dir/../run_output
[ -d $output_dir ] || mkdir $output_dir

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results 
[ -d $result_dir ] || mkdir $result_dir


########################
## 1. QC/mapping/assembly/quantification for all samples in the input dir  (Tophat/Cufflink/Htseq-count)
########################
cd $input_dir

c=0;h=0;gtflist="";samlist=""; labels=""

for i in *.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q big-multi $cpu $memory $email _RNAseq.sh $R1 $R2;

done

exit

########################
## 2. merge all samples to get big matrix for expression (e.g. one row per gene/Tx, one col per sample)
########################
[ -d $result_dir/merged ] || mkdir $result_dir/merged
cd $result_dir/merged

Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/uniq/genes.fpkm_tracking` genes.fpkm.allSamples.tab
Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/uniq/isoforms.fpkm_tracking` isoforms.fpkm.allSamples.tab

########################
## 3. draw aggregation plot for the RNAseq density in the genetic region
# Note: see note in bigWigAverageOverBed_81bins.sh for requisition 
########################
[ -d $result_dir/aggregationPlot ] || mkdir $result_dir/aggregationPlot
cd $result_dir/aggregationPlot
LOG=$result_dir/aggregationPlot/_aggPlot.log
for i in $fordisplay_dir/*uniq*ed.bw;
do
    bsub -J aggPlot -oo $LOG.$i -eo $LOG.$i -q big-multi $email -M 4000 -R rusage[mem=4000] -n 2 aggregationPlot_getBinSignal.sh $GENOME/Annotation/Genes/gencode.v19.annotation.bed12.mRNA.81bins.bed 81 $i;
done

# set break point here to wait until all samples are completedly processed.
Rscript $pipeline_path/src/draw_aggregation_plot.R

########################
## 4. Added cluster procedure -- by Bin
########################
bsub Rscript _clustComRNASeq.R $output_dir $result_dir

###########
# RNA-SeQC
#
# make samplelist_file.txt

./makeSamplelistFile.sh "$output_dir/*/acce*.bam" > $result_dir/RNA-SEQCfolder/samplelist_file.txt

# gencodePlusMask.v13.annotation.gtf is a combination of 
# gencode.v13.annotation.karotyped.gtf and chrM.rRNA.tRNA.gtf

bsub -J runRnaSeqc -q big-multi -n 4 -R 'rusage[mem=10000]' "java -jar -Xmx64g RNA-SeQC_v1.1.7.jar -s $result_dir/RNA-SEQCfolder/samplelist_file.txt -t $result_dir/RNA-SEQCfolder/gencodePlusMask.v13.annotation.gtf -r $GENOME/Sequence/Bowtie2Index/genome.fa -o $result_dir/RNA-SEQCfolder/RNASeQCoutput"


########################
## 5. factor analysis to identify the hidden covariates (PEER)
########################
bsub Rscript _factor_analysis.R

########################
## 6. identify differentially expressed genes (cuffdiff and DEseq), incoperating the hidden covariates from PEER
########################
[ -d $result_dir/DE_cuffdiff ] || mkdir $result_dir/DE_cuffdiff
cd $result_dir/DE_cuffdiff

bsub -o _DE_cuffdiff.log -q long $cpu $memory $email _DE_cuffdiff.sh $gtflist $samlist $labels

[ -d $result_dir/DE_DESeq2 ] || mkdir $result_dir/DE_DESeq2
cd $result_dir/DE_DESeq2
# Save the covariance table from Google Doc
# wget --no-check-certificate -qO - "https://docs.google.com/spreadsheet/ccc?key=0Aumm3V3g3dF7dEFnZ2pPQjlheXlZand6YWUxeF9PMUE&gid=5&output=txt" > covariances.tab
bsub -o _DE_DESeq2.log -q long $cpu $memory $email Rscript _DE_DESeq2.R $output_dir PD Ct $ANNOTATION

########################
## 7. eQTL 
########################
[ -d $result_dir/eQTL ] || mkdir $result_dir/eQTL
cd $result_dir/eQTL

## prepare input files
## ===================
# SNP (call from the unique mapper)
for i in $output_dir/*/uniq/accepted_hits.snp.depth_gt_15; do echo $i; cut -f1-4 $i > $i.bdg; done
bedtools unionbedg -i $output_dir/*/uniq/accepted_hits.snp.depth_gt_15.bdg -filler -:- > union.snp  # 299003 in total
awk 'gsub("-:-","-:-") <= 0.95*(NF-3)' union.snp > union.MAF5p.snp  # 20698 in total with MAF>=0.05
#awk 'gsub("-:-","-:-") <= 0.5*(NF-3)' union.snp > union.common38.snp  # 2413 in total occurring in 50% samples
# get reads depth per allele
samtools mpileup -Bf $GENOME/Sequence/Bowtie2Index/genome.fa -l union.MAF5p.snp $output_dir/*/uniq/accepted_hits.bam > union.pileup  
# parse the pileup file to get allele depth and then encode into additive code
cat union.pileup | awk --posix -f $pipeline_path/modules/_pileup2depth.awk | awk -f $pipeline_path/modules/_depth2additive.awk > tmp1

# annotate SNP
dbSNP=$GENOME/Annotation/Variation/snp137.bed.groupped.SNP
awk '{OFS="\t";split($1, s, "_"); print s[1],s[2]-1, s[2], $1;}' tmp1 | bedtools sort | intersectBed -a - -b $dbSNP  -sorted -wao | cut -f1-4,8 | groupBy -g 1,2,3,4 -c 5 -o collapse | sort -k4,4 | cut -f4-5 | paste - <(sort -k1,1 tmp1) | awk '{OFS="\t"; $3=$3"_"$2; print}' | cut -f3- > tmp2

echo -ne "sampleName\t" > tmp1
for i in $output_dir/*/uniq/accepted_hits.bam; do echo $i | sed 's/.*run_output\///g;s/\/uniq.*//g';done | rowsToCols stdin stdout >> tmp1
# note: some snp might have MAF<0.05 even it has >16 reads supported
awk '{s=0;for(i=2;i<=NF;i++) {if($i>0) s++;} if(s/(NF-1) >= 0.05) print}' tmp2 | cat tmp1 - > snp.additive.txt

# exclude chrM and only *SNDA_[2-4] samples
grep -v chrM snp.additive.txt | rowsToCols stdin stdout | grep -E "sample|SNDA_[2-4]" | rowsToCols stdin snp.txt

# expression
#cat $result_dir/DE_DESeq2/PDvsHC/htseqcount.vst.allsamples.xls | rowsToCols stdin stdout | grep -E "sample|SNDA_[2-4]" | rowsToCols stdin expression.txt
# 
cat $result_dir/merged/genes.fpkm.allSamples.tab | rowsToCols stdin stdout | grep -E "tracking_id|gene_short_name|SNDA_[2-4]" | rowsToCols stdin stdout | sed 's/\t/__/' | awk '{if(NR==1) {print; next;} s=0;for(i=2;i<=NF;i++)s=+$i; if(s/(NF-1) > 0.01) print}' > expression.txt

# covariance
grep -E "sample|SNDA_[2-4]" $input_dir/covariances.tab | rowsToCols stdin stdout | sed 's/\bHC\b/0/g;s/\bILB\b/1/g;s/\bPD\b/2/g;s/\bF\b/0/g;s/\bM\b/1/g;s/\bSNDA\b/0/g;s/\bMCPY\b/1/g;s/\bTCPY\b/2/g' | grep -v "cellType" > cov.txt

# gene position
echo -e "geneid\tchr\ts1\ts2" > geneloc.txt
#cut -f1,7 $result_dir/merged/genes.fpkm.allSamples.tab | tail -n +2 | sed 's/[:-]/\t/g' >> geneloc.txt
cut -f1,5,7 $result_dir/merged/genes.fpkm.allSamples.tab | tail -n +2 | sed 's/\t/__/' | awk '{OFS="\t"; gsub("[:-]","\t", $2); print}' >> geneloc.txt
# SNP position
echo -e "snp\tchr\tpos" > snpsloc.txt
cut -f1 snp.txt | tail -n +2 | awk '{OFS="\t";split($1, s, "_"); print $1,s[1],s[2];}' >> snpsloc.txt

# run eQTL
#Rscript $pipeline_path/modules/_eQTL.R snp.txt expression.txt cov.txt output
Rscript $pipeline_path/modules/_eQTL.R snp.txt expression.txt cov.txt output geneloc.txt snpsloc.txt
## plot expression vs. genotype
Rscript $pipeline_path/modules/_eQTL.plot.R output.cis.txt snp.txt expression.txt cov.txt

grep -vE "cellType|condition" cov.txt > cov2.txt
Rscript $pipeline_path/modules/_eQTL.R snp.txt expression.txt cov2.txt output_wtcondition geneloc.txt snpsloc.txt
Rscript $pipeline_path/modules/_eQTL.plot.R output_wtcondition.cis.txt snp.txt expression.txt cov.txt
########################
## 8. pathway analysis (SPIA)
########################
Rscript _pathway_analysis.R
