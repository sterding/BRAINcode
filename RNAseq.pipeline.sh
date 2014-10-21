####################################
# Pipeline for RNAseq data Analysis
# Authors: Bioinformatics Team @ Scherzer's lab 
# Email: xdong@rics.bwh.harvard.edu
# Date: 10/6/2014
# Version: 2.0
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
pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

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
## 1. Processing per sample 
########################
cd $input_dir

c=0;h=0;gtflist="";samlist=""; labels=""

for i in *.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}
    
    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] -u $EMAIL -N _RNAseq.sh $R1 $R2;
    #bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q $QUEUE -n $CPU -M $MEMORY -u $EMAIL -N _RNAseq.sh $R1 $R2;

done

exit

########################
## [test] (merge all reads and then call assembly once) vs. (call assembly individually and then run cuffmerge)
########################
samtools merge -1r all124samples.uniq.accepted_hits.bam `ls */uniq/accepted_hits.bam.non-rRNA-mt.bam`


########################
## [test] total-rRNA-chrM vs. total, which one is better to be used for normalization?
########################

echo "sampleID" `rowsToCols /PHShome/xd010/neurogen/rnaseq_PD/run_output/PD_UWA734_SNDA_2/uniq/accepted_hits.bam.bam2annotation stdout | sed 's/://g' | head -n1` | sed 's/ /\t/g' > allsamples.bam2annotation.tab
paste <(ls -1 ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.bam2* | sed 's/.*run_output\///g;s/\/uniq.*//g') <(paste ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.bam2* | sed 's/[^[:space:]]*: //g' | rowsToCols stdin stdout) >> allsamples.bam2annotation.tab
#R
pdf("allsamples.bam2annotation.pdf")
df=read.table("allsamples.bam2annotation.tab", header=T)
rownames(df)=df[,1]; df=df[,-1]; df=df/1e6
plot(total_non.rRNA.mt~total, df, asp=1, xlim=c(0,300), ylim=c(0,300), col=gsub(".*_([1-5]).*","\\1",rownames(df)), xlab="total reads (in million)", ylab="total-rRNA-chrM (in million)")
abline(a=0,b=1,lty=2)
legend("topleft", paste("batch", 1:5), col=1:5, bty='n', pch=1)
dev.off()

########################
## 2. merge all samples to get big matrix for expression (e.g. one row per gene/Tx, one col per sample)
########################
[ -d $result_dir/merged ] || mkdir $result_dir/merged
cd $result_dir/merged

# multi mapper
Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/genes.fpkm_tracking` genes.fpkm.allSamples.multi.xls
Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/isoforms.fpkm_tracking` isoforms.fpkm.allSamples.multi.xls
Rscript $pipeline_path/modules/_mergeSamples_htseq.R `ls $output_dir/*/hgseqcount.by.gene.tab` genes.htseqcount.allSamples.multi.xls

# uniq mapper
Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/uniq/rpkm/genes.fpkm_tracking` genes.fpkm.allSamples.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples.R `ls $output_dir/*/uniq/rpkm/isoforms.fpkm_tracking` isoforms.fpkm.allSamples.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples_htseq.R `ls $output_dir/*/uniq/hgseqcount.by.gene.tab` genes.htseqcount.allSamples.uniq.xls

#--------------------------
# 2.2 combined bigwig into 
#--------------------------
for i in HC_TCPY HC_MCPY HC_SNDA ILB_SNDA PD_SNDA;
do
    [ "$i" = "HC_TCPY" ]  && pattern="(HC|ND)_.*_TCPY_[2345]";
    [ "$i" = "HC_MCPY" ]  && pattern="(HC|ND)_.*_MCPY_[2345]";
    [ "$i" = "HC_SNDA" ]  && pattern="(HC|ND)_.*_SNDA_[2345]";
    [ "$i" = "ILB_SNDA" ] && pattern="ILB_.*_SNDA_[2345]";
    [ "$i" = "PD_SNDA" ]  && pattern="PD_.*_SNDA_[2345]";
    
    echo $i, $pattern;
    
    bsub -J combine_bw -oo _combin_bw.log -eo _combin_bw.log -q $QUEUE -n $CPU -M $MEMORY -u $EMAIL -N _combine_bigwig.sh $i $pattern
done

#--------------------------
## 2.3 make UCSC track hub
#--------------------------
_make_trackDb.sh > $fordisplay_dir/trackDb.RNAseq.txt
rsync -azvL $fordisplay_dir/trackDb.RNAseq.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
rsync -azvL $fordisplay_dir/*uniq.accepted_hits.normalized2.bw xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/
rsync -azv *.xls *.bw xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/merged

#--------------------------
### expression trend of specific gene(s) along the stages (e.g. HC,ILB,PD)
#--------------------------
grep -w -P "tracking_id|SNCA" isoforms.fpkm.allSamples.uniq.xls | Rscript $pipeline_path/modules/_plotTrend.R stdin SNCA.tx.pdf
grep -w -P "tracking_id|SNCA|GBA|LRRK2" genes.fpkm.allSamples.uniq.xls | Rscript $pipeline_path/modules/_plotTrend.R stdin SNCA.pdf

#########################
### 3. use RLE to detect outlier
#########################
Rscript $pipeline_path/modules/_normQC.R genes.fpkm.allSamples.uniq.xls genes.fpkm.allSamples.uniq.QC.pdf

########################
## 4. factor analysis to identify the hidden covariates (PEER)
########################
Rscript _factor_analysis.R

#########################
### 5. post-normalization QC (use RLE to detect outlier)
#########################
## TODO:
# 

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
