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
  echo "Usage:
  $HOME/neurogen/pipeline/RNAseq/RNAseq.pipeline.sh /data/neurogen/rnaseq_PD/rawfiles
  "
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


## TODO: test the prerequisitions, incl.
# kpal, fastqc, tophat, bowtie, CIRI, GATK, cufflinks, htseq-count, bedtools, samtools, R, fastq-mcf, python
# Jim-Kent's utility: bigWigSummary ...
# R package: DESeq2, MatrixEQTL, SPIA, SVA, PEER, ggplot2 etc.

########################
## 1. Processing per sample
########################
cd $input_dir

for i in *.R1.fastq.gz;
do
    R1=$i
    R2=${i/R1/R2};
    samplename=${R1/.R1*/}

    # run the QC/mapping/assembly/quantification for RNAseq
    bsub -J $samplename -oo $output_dir/$samplename/_RNAseq.log -eo $output_dir/$samplename/_RNAseq.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] -u $EMAIL -N _RNAseq.sh $R1 $R2;

done

exit

########################
## [test] total-rRNA-chrM vs. total, which one is better to be used for normalization?
########################
[ -d $result_dir/coverage ] || mkdir $result_dir/coverage
cd $result_dir/coverage

echo "sampleID" `rowsToCols /PHShome/xd010/neurogen/rnaseq_PD/run_output/PD_UWA734_SNDA_2_rep1/uniq/accepted_hits.bam.bam2annotation stdout | sed 's/://g' | head -n1` | sed 's/ /\t/g' > allsamples.bam2annotation.tab
paste <(ls -1 ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.bam2annotation | sed 's/.*run_output\///g;s/\/uniq.*//g') <(paste ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.bam2annotation | sed 's/[^[:space:]]*: //g' | rowsToCols stdin stdout) >> allsamples.bam2annotation.tab
#R
pdf("allsamples.bam2annotation.pdf")
df=read.table("allsamples.bam2annotation.tab", header=T)
rownames(df)=df[,1]; df=df[,-1]; df=df/1e6
plot(total_non_rRNA_mt~total, df, asp=1, xlim=c(0,300), ylim=c(0,300), col=gsub(".*_([1-9])_.*","\\1",rownames(df)), xlab="total reads (in million)", ylab="total-rRNA-chrM (in million)")
abline(a=0,b=1,lty=2)
legend("topleft", paste("batch", sort(unique(as.numeric(gsub(".*_([1-9])_.*","\\1",rownames(df)))))), col=sort(unique(as.numeric(gsub(".*_([1-9])_.*","\\1",rownames(df))))), bty='n', pch=1)

## group them into HC_TCPY HC_MCPY HC_SNDA ILB_SNDA PD_SNDA HC_SNDA HCILB_SNDA HC_PBMC HC_FB HC_SN HC_SNDAstranded

dev.off()

########################
# 2.2 calculatinng coverage of RNAseq
########################
[ -d $result_dir/coverage ] || mkdir $result_dir/coverage
cd $result_dir/coverage

## cumulative coverage for all samples in the group
# ----------------------------------
for i in HCILB_SNDA HC_PY HC_nonNeuron;
do
    #bsub -J combine_coverage -oo _combine_cov.$i.log -eo _combine_cov.$i.log -q $QUEUE -n $CPU -M $MEMORY -u $EMAIL -N _combine_coverage.sh $i 5
    bsub -q big-multi -n 4 -M 8000 _combine_coverage.sh $i 0.05
done

## cumulative coverage for 7 samples in each group
# ----------------------------------
for i in `seq 1 100`;
do
  bsub -q normal -n 2 -M 2000 _combine_coverage.sh HCILB_SNDA 0.05 7; # output: random.covered.0.05RPM.HCILB_SNDA.*.txt
  bsub -q normal -n 2 -M 2000 _combine_coverage.sh HC_PY 0.05 7;
done
cut -f2 random.covered.0.05RPM.HCILB_SNDA.*.txt | paste - - - - - - - > covered.0.05RPM.HCILB_SNDA.random7.txt
cut -f2 random.covered.0.05RPM.HC_PY.*.txt | paste - - - - - - - > covered.0.05RPM.HC_PY.random7.txt
rm random.covered.0.05RPM.HC*

##  coverage for individual sample with different RPM cutoff (>0, >=0.01, >=0.05, >=0.1, >=0.5, >=1)
# ----------------------------------
for i in ../../run_output/*/uniq/accepted_hits.normalized.bedGraph; do echo $i; bsub -q short -n 1 -M 500 "awk 'BEGIN{s0=0;s1=0;s2=0;s3=0;s4=0;s5=0;}{L=(\$3-\$2); (\$4>=1)?(s5+=L):((\$4>=0.5)?(s4+=L):((\$4>=0.1)?(s3+=L):((\$4>=0.05)?(s2+=L):((\$4>=0.01)?(s1+=L):(s0+=L)))));}END{print s0,s1,s2,s3,s4,s5;}' $i > $i.coverageWithRPM"; done
echo -e "sample\tRPMgt0\tRPMgt0.01\tRPMgt0.05\tRPMgt0.1\tRPMgt0.5\tRPMgt1" > coverageWithRPM.txt
paste <(ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph.coverageWithRPM | sed 's/.*run_output\/\(.*\)\/uniq.*/\1/g') <(cat ../../run_output/*/uniq/accepted_hits.normalized.bedGraph.coverageWithRPM) | sed 's/ /\t/g' >> coverageWithRPM.txt

## body coverage
module load RSeQC/2.4
for i in ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.non-rRNA-mt.bam; do bsub -J $i -q short -n 1 -M 3000 -R 'rusage[mem=3000]' geneBody_coverage.py -r $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 -i $i  -o $i; done
cat `ls ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.non-rRNA-mt.bam.geneBodyCoverage.txt | grep -f ../merged/samplelist.HCILB_SNDA` | grep accepted_hits > geneBodyCoverage.HCILB_SNDA.txt

## cumulative coverage
for i in ~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA_*rep[1-2]/uniq/accepted_hits.bw;
do
    echo $i;
    [ -f $i.gt5reads.bedGraph ] || bigWigToBedGraph $i stdout | awk '$4>=5' > $i.gt5reads.bedGraph &
done
bsub -J combine_coverage -oo _combine_cov.HCILB_SNDA.log -eo _combine_cov.HCILB_SNDA.log -q $QUEUE -n $CPU -M $MEMORY -u $EMAIL -N _combine_coverage.sh HCILB_SNDA 5

> tmp.bed
> cummulative.cov.txt
for i in `ls -1 ~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA_*rep[1-2]/uniq/accepted_hits.bw.gt5reads.bed -S`;
do
    awk -vi=$i '{s+=($3-$2)}END{print i, s, s/3137161264}' $i;
    cat tmp.bed $i | sortBed | mergeBed -i - > tmp2.bed
    awk -vi=$i '{s+=($3-$2)}END{print i, s, s/3137161264}' tmp2.bed >> cummulative.cov.txt
    mv tmp2.bed tmp.bed
done

## Rscript to make all kinds of coverage plots
Rscript $HOME/neurogen/pipeline/RNAseq/module/_coverage.plots.R

## reads count for
echo 'ID' `sort ~/neurogen/rnaseq_PD/run_output/PD_UWA734_SNDA_6_rep2/uniq/accepted_hits.bam.bam2annotation | cut -f1 -d':' | rowsToCols stdin stdout` | sed 's/ /\t/g' > uniq.bam2annotation.stat.txt
for i in ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam.bam2annotation; do echo $i `sort $i | cut -f2 -d' ' | rowsToCols stdin stdout`; done | sed 's/.*run_output\/\(.*\)\/uniq.*tion/\1/g;s/ /\t/g' >> uniq.bam2annotation.stat.txt
# import uniq.bam2annotation.stat.txt into Google spreadsheet: https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/edit#gid=289128536

########################
## 2. merge all samples to get big matrix for expression (e.g. one row per gene/Tx, one col per sample)
########################
[ -d $result_dir/merged ] || mkdir $result_dir/merged
cd $result_dir/merged

# uniq mapper
Rscript $pipeline_path/modules/_mergeSamples.R `ls ../../run_output/*/uniq/rpkm/genes.fpkm_tracking` genes.fpkm.cufflinks.allSamples.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples.R `ls ../../run_output/*/uniq/rpkm/isoforms.fpkm_tracking` isoforms.fpkm.cufflinks.allSamples.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples_htseq.R `ls ../../run_output/*/uniq/hgseqcount.by.gene.tab` genes.htseqcount.cufflinks.allSamples.uniq.xls

# HC/ILB only
Rscript $pipeline_path/modules/_mergeSamples.R `ls ../../run_output/*/uniq/rpkm/genes.fpkm_tracking | grep -v PD_` genes.fpkm.cufflinks.HCILB.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples.R `ls ../../run_output/*/uniq/rpkm/isoforms.fpkm_tracking | grep -v PD_` isoforms.fpkm.cufflinks.HCILB.uniq.xls
Rscript $pipeline_path/modules/_mergeSamples_htseq.R `ls ../../run_output/*/uniq/hgseqcount.by.gene.tab | grep -v PD_` genes.htseqcount.cufflinks.HCILB.uniq.xls

# intronic FPKM
cut -f1 ../../run_output/PD_UWA734_SNDA_6_rep2/uniq/meanRPM.of.metaintron.by.gene.tab | grep ENSG | sort > metaIntron.meanRPM.allSamples.xls
for i in ../../run_output/*/uniq/meanRPM.of.metaintron.by.gene.tab;
do
  echo $i;
  grep ENSG $i | sort | cut -f4 | paste metaIntron.meanRPM.allSamples.xls - > /tmp/metaIntron.meanRPM.allSamples.xls
  cp /tmp/metaIntron.meanRPM.allSamples.xls metaIntron.meanRPM.allSamples.xls
done
echo "locus" `ls ../../run_output/*/uniq/meanRPM.of.metaintron.by.gene.tab | sed 's/.*run_output\/\(.*\)\/uniq.*/\1/g' | rowsToCols stdin stdout` | sed 's/ /\t/g' | cat - /tmp/metaIntron.meanRPM.allSamples.xls > metaIntron.meanRPM.allSamples.xls

## UPdate: use cuffquant --> cuffnorm to calculate normalized expression FPKM
bsub -J cuffnorm -oo _cuffnorm.log -eo _cuffnorm.log -q big-multi -n 8 -M 10000 -R rusage[mem=10000] cuffnorm -o ./cuffnorm --no-update-check -L `ls /data/neurogen/rnaseq_PD/run_output/*/uniq/rpkm/abundances.cxb | sed 's/.*run_output\/\(.*\)\/uniq.*/\1/g' | tr '\n' ','` -p 8 -total-hits-norm -library-norm-method quartile $ANNOTATION_GTF `ls /data/neurogen/rnaseq_PD/run_output/*/uniq/rpkm/abundances.cxb`
for i in cuffnorm/*.count_table cuffnorm/*.fpkm_table; do echo $i; sed 's/_0//g' $i > $i.tmp; mv $i.tmp $i; done
ln -fs cuffnorm/genes.fpkm_table genes.fpkm.cuffnorm.allSamples.uniq.xls
ln -fs cuffnorm/genes.count_table genes.count.cuffnorm.allSamples.uniq.xls
ln -fs cuffnorm/isoforms.fpkm_table isoforms.fpkm.cuffnorm.allSamples.uniq.xls
ln -fs cuffnorm/isoforms.count_table isoforms.count.cuffnorm.allSamples.uniq.xls

#--------------------------
# 2.2 combined bigwig into
#--------------------------
for i in HCILB_SNDA HC_PY HC_nonNeuron HC_Neuron HC_TCPY HC_MCPY HC_SNDA ILB_SNDA PD_SNDA HC_SN HC_PBMC HC_FB;
do
    [ -e trimmedmean.uniq.normalized.$i.bw ] || bsub -J combine_bw -oo _combin_bw.$i.log -eo _combin_bw.$i.log -q normal -n 4 -M 6000 -R rusage[mem=6000] _combine_bigwig.sh $i
done


#--------------------------
## 2.3 make UCSC track hub
#--------------------------
_make_trackDb.sh > $fordisplay_dir/trackDb.RNAseq.txt
rsync -azvL $fordisplay_dir/trackDb.RNAseq.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
rsync -azvL $fordisplay_dir/*uniq.accepted_hits.normalized.bw xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/
rsync -azv *.xls *.bw xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/merged

#--------------------------
### expression trend of specific gene(s) along the stages (e.g. HC,ILB,PD)
#--------------------------
grep -w -P "tracking_id|SNCA" isoforms.fpkm.cufflinks.allSamples.uniq.xls | Rscript $pipeline_path/modules/_plotTrend.R stdin SNCA.tx.pdf
grep -w -P "tracking_id|SNCA|GBA|LRRK2" genes.fpkm.cufflinks.allSamples.uniq.xls | Rscript $pipeline_path/modules/_plotTrend.R stdin SNCA.pdf


#########################
### 3. detect outlier
#########################
cd $result_dir/merged
## k-mer distance
kpal cat $output_dir/[!PD]*/k9 mergedHCILB_k9 # not the controls (e.g. stranded, amplified)
kpal matrix mergedHCILB_k9 mergedHCILB_k9.matrix.txt

Rscript $pipeline_path/modules/_normQC.R genes.fpkm.cufflinks.HCILB.uniq.xls mergedHCILB_k9.matrix.txt QC.genes.fpkm.cufflinks.HCILB.uniq.pdf
Rscript $pipeline_path/modules/_normQC.R genes.fpkm.cuffnorm.HCILB.uniq.xls mergedHCILB_k9.matrix.txt QC.genes.fpkm.cuffnorm.HCILB.uniq.pdf

#########################
### 3. check consistency of replicates
#########################
cd $result_dir/merged
Rscript $pipeline_path/modules/_replicates.R genes.fpkm.cufflinks.HCILB.uniq.xls QCrep.genes.fpkm.HCILB.uniq.pdf

# linear amplification vs. non-amplification
#Rscript $pipeline_path/modules/_pairwise_compare.R HC_M0235-4_PBMC_6_rep1.amplified HC_M0235-4_PBMC_6_rep1.unamplified
Rscript $pipeline_path/modules/_pairwise_compare.R HC_UWA616_SN_6_rep1.amplified HC_UWA616_SN_6_rep1.unamplified

# laser-captured neurons vs. homogeneous brain tissue (see the code at the end of _pairwise_compare.R)
#Rscript $pipeline_path/modules/_pairwise_compare.R HC_UWA616_SN_6_rep1.amplified HC_UWA616_SNDA_2_rep1 

# same sample, same batch
Rscript $pipeline_path/modules/_pairwise_compare.R HC_BN05-10_SNDA_5_rep1 HC_BN05-10_SNDA_5_rep2
#Rscript $pipeline_path/modules/_pairwise_compare.R HC_BN08-90_SNDA_5_rep1 HC_BN08-90_SNDA_5_rep2

# same sample, different batch
#Rscript $pipeline_path/modules/_pairwise_compare.R HC_MGH1000_SNDA_1_rep1 HC_MGH1000_SNDA_4_rep2
# Rscript $pipeline_path/modules/_pairwise_compare.R PD_MGH1288_SNDA_1_rep1 PD_MGH1288_SNDA_4_rep2
# Rscript $pipeline_path/modules/_pairwise_compare.R PD_MGH1488_SNDA_1_rep1 PD_MGH1488_SNDA_4_rep2
# Rscript $pipeline_path/modules/_pairwise_compare.R PD_UWA734_SNDA_2_rep1 PD_UWA734_SNDA_6_rep2

# brain vs. non-brain
# SN vs. PBMC
Rscript $pipeline_path/modules/_pairwise_compare.R HC_UWA616_SN_6_rep1.amplified HC_B0254-4_PBMC_6_rep1
# SN vs. FB
Rscript $pipeline_path/modules/_pairwise_compare.R HC_UWA616_SN_6_rep1.amplified HC_ND34770_FB_6_rep1


########################
## 4. eQTL
########################
cd $result_dir/eQTL

# ## eQTL with PEER normalization (DEPRECATED)
# module unload R/3.1.0
# module load R/3.0.2
# Rscript $pipeline_path/modules/_PEER_eQTL.R $result_dir/merged/genes.fpkm.cufflinks.HCILB.uniq.xls 
# module unload R/3.0.2; module load R/3.1.0


## eQTL with SVA normalization
cd $result_dir/eQTL/HCILB_SNDA
bsub -q big -n 2 -R 'rusage[mem=10000]' -eo eQTL.run.log -oo eQTL.run.log Rscript ~/neurogen/pipeline/RNAseq/modules/_SVA.eQTL.R  
ln -fs final.cis.eQTL.xls final.cis.eQTL.d1e6.p1e-2.xls
cat final.cis.eQTL.xls | awk '$6<=0.05' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls

## run permutations in bash
mkdir permutations
for i in `seq 1 10000`; do [ -e permutations/permutation$i.txt ] || bsub -n 1 -M 500 -q short -J $i Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_permutation_minP.R $i data.RData expression.postSVA.xls; done

# post-eQTL analysis
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.R

## RTC analysis for eQTL (output: final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.xls with 11 columns)
bsub -q big-multi -n 4 -M 10000 -oo RTC.run.log -eo RTC.run.log Rscript $pipeline_path/modules/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls gene
## annotate: add hostgene_GWAS_SNP and hostgene_eQTL_SNP
sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.xls | awk '{OFS="\t"; split($9,a,"_"); if(NR>1) print a[1],$10-1,$10,$0}'  | intersectBed -a - -b <(cat $GENOME/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-14,18 | sort | groupBy -g 1-11 -c 12 -o distinct | awk '{OFS="\t"; split($9,a,"_"); print a[1],$11-1,$11,$0}' | intersectBed -a - -b <(cat $GENOME/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-15,19 | sort | groupBy -g 1-12 -c 13 -o distinct | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.annotated.xls
## filter: for each gene, take the eSNP with the best RTC score (if there are multiple eSNPs in LD) per GWAS SNP
sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.annotated.xls | sort -k2,2r -k5,5 -k7,7gr | awk '{if(id!=$2$5) print; id=$2$5;}' | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls

## manhanten plot for eQTL with RTC
# Note: Run RTC (_RTC.R) ahead to get final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls
Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls mRNA
Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls ncRNA

# simplfied table with LD infor (take the eSNP with the smallest p-value per SNP/associatedGene/Trait)
# Note: we use LD block called by "PLINK --blocks", which use hyplotypeviewer behind and different from pairwise LD caller like SNAP)
echo "#SNP OmniID Ref:Alt Chr eSNP_host_gene Associated_transcript_hostgene Associated_transcript minP Trait RTC GWAS_SNP GWAS_SNP_pos GWAS_SNP_pvalue LD" | sed 's/ /\t/g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.simplified.xls
sed 's/ /_/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls | awk '{OFS="\t"; if($7>=0.85) {split($9,a,"_");print a[1],$11-1,$11,$1,a[1],$13,($8=="NA")?"(intergenic)":$8,$9,$3,$6,$7,$5,$10}}' | intersectBed -a - -b ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed -wo | cut -f4-13,17 | awk '{OFS="\t"; $1=$1"|"$11; print}' | cut -f1-10 | sort -k1,1 -k3,3 -k4,4 -k7,7 -k6,6g | awk '{if(id!=$1$3$4$7) {print; id=$1$3$4$7;}}' | while read rs chr rest; do chr=${chr/chr/}; rs0=${rs/\|*/}; ld=`fgrep -w $rs0 $GENOME/Annotation/Variation/1000G/LDblock/Chr$chr.LD.blocks.det | head -n1 | awk '{print $1"_"$2"_"$3}'`; echo $rs $chr $rest chr$ld; done | sed 's/ /\t/g' | awk '{OFS="\t"; print "chr"$2,$10-1,$10,$9,$0}' | intersectBed -a - -b <(sed 's/ /_/g' /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pvalue.bed) -wo | awk '$11==$22' | cut -f5-15,20 | awk '{OFS="\t"; split($1,a,"[|_]"); print a[1],a[2],a[3],$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$11}' >> final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.simplified.xls

## Top eQTL at 1M bp stepping
grep -v protein_coding $GENOME/Annotation/Genes/genes.bed | cut -f4 | fgrep -w -f - ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls | awk '{OFS="\t"; if($5<=1e-6) print $8,$7-1,$7,$1"|"$2,-log($5)/log(10),-log($6)/log(10);}' | sortBed | intersectBed -a - -b $GENOME/Annotation/Genes/genes.bed -wo | cut -f1-2,4-6,13-14 | awk '{OFS="\t"; print int($2/1000000),$0}' | sort -k2,2 -k1,1n -k5,5gr | awk '{OFS="\t"; if(id!=$1$2) {id=$1$2;p=$5;print;} else if($5==p) print;}' | sed 's/|/\t/g' | sort -k1,1 -k2,2 -k5,5 -k6,6 -k8,8 | groupBy -g 1,2,5-9 -c 3,3,4 -o count,distinct,distinct -delim "|" | sort -k3,3 | join -1 3 -2 4 - <(sort -k4,4 $GENOME/Annotation/Genes/genes.bed) -a 1 | awk '{OFS="\t"; print "ncRNA",$2,$3,$1,$16"___"$1"___"$17,$4,$5,$6,$7,$8,$9,$10}' | sort -k3,3 -k2,2n -k6,6gr > ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt
grep protein_coding $GENOME/Annotation/Genes/genes.bed | cut -f4 | fgrep -w -f - ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls | awk '{OFS="\t"; if($5<=1e-6) print $8,$7-1,$7,$1"|"$2,-log($5)/log(10),-log($6)/log(10);}' | sortBed | intersectBed -a - -b $GENOME/Annotation/Genes/genes.bed -wo | cut -f1-2,4-6,13-14 | awk '{OFS="\t"; print int($2/1000000),$0}' | sort -k2,2 -k1,1n -k5,5gr | awk '{OFS="\t"; if(id!=$1$2) {id=$1$2;p=$5;print;} else if($5==p) print;}' | sed 's/|/\t/g' | sort -k1,1 -k2,2 -k5,5 -k6,6 -k8,8 | groupBy -g 1,2,5-9 -c 3,3,4 -o count,distinct,distinct -delim "|" | sort -k3,3 | join -1 3 -2 4 - <(sort -k4,4 $GENOME/Annotation/Genes/genes.bed) -a 1 | awk '{OFS="\t"; print "mRNA",$2,$3,$1,$16"___"$1"___"$17,$4,$5,$6,$7,$8,$9,$10}' | sort -k3,3 -k2,2n -k6,6gr >> ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt

## manhanten plot for eQTL
Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls mRNA
Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls ncRNA

## boxplot of top eQTL
cut -f4,12 ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt | cut -f1 -d"|" > topeQTL.gene.snp.list
Rscript ~/neurogen/pipeline/RNAseq/src/eQTLlist2plot.R topeQTL.gene.snp.list

## boxplot of top RTC eQTL

cat final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls | awk '{FS="\t"; OFS="\t"; split($9,a,"_"); if(NR>1) print a[1],$11-1,$11,$1,$2}' | sortBed | intersectBed -a - -b ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed -wo| cut -f5,9 > RTC.gene.snp.list
Rscript ~/neurogen/pipeline/RNAseq/src/eQTLlist2plot.R RTC.gene.snp.list


## disease SNP enrichment in the eQTL SNPs
## =====================================
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
SNPpos=~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID
awk '$6<=0.05' final.cis.eQTL.xls | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.txt
awk '$6<=0.05' final.cis.eQTL.xls | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.txt

# split into mRNA and ncRNA
awk '$6<=0.05' final.cis.eQTL.xls | fgrep -wf <(fgrep -w protein_coding /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | cut -f4) | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.mRNA.txt
awk '$6<=0.05' final.cis.eQTL.xls | fgrep -wf <(fgrep -v protein_coding /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | cut -f4) | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.ncRNA.txt
awk '$6<=0.05' final.cis.eQTL.xls | fgrep -wf <(fgrep -w protein_coding /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | cut -f4) | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.mRNA.txt
awk '$6<=0.05' final.cis.eQTL.xls | fgrep -wf <(fgrep -v protein_coding /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed | cut -f4) | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.ncRNA.txt


R
setwd("~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA")
N=as.numeric(system("wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '", intern=T)) ## total SNPs in dbSNP: 48709140
nDiseases=as.numeric(system("cut -f4 $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed.autosomal.associations.bed | sort -u | wc -l | cut -f1 -d' '", intern=T)) ## total diseasese in gwas
n=as.numeric(system("more +2 final.cis.eQTL.xls | awk '$6<=0.05' | cut -f1,7 | sort -u | wc -l", intern=T))  ## total eSNP
results=data.frame();
for(i in c('HTNE','gene','mRNA','ncRNA')){
  df1=read.table(paste0("dSNP.with.eSNP.",i,".txt"), header=F); rownames(df1)=df1[,1]
  df2=read.table(paste0("dSNP.without.eSNP.",i,".txt"), header=F); rownames(df2)=df2[,1]
  common=intersect(rownames(df1), rownames(df2))
  df=cbind(df1[common,], df2[common,2]); df=df[,-1]; colnames(df)=c('dSNP_eSNP','dSNP_N_eSNP')  ## only the disease with both eSNP and dSNP
  results = rbind(results,
                  cbind(Disease_or_Trait=rownames(df),
                  df,
                  type=i,
                  pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$p.value),
                  OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$estimate))
                )
}
results=subset(results, OR>1 & pvalue<0.05/nDiseases & dSNP_eSNP>3)
results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
results = results[with(results, order(type, -dSNP_eSNP)), ]

write.table(results, "dSNP.eQTL.enrichment.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)

results$pvalue[results$pvalue==0]=2.2e-16
results = results[with(results, order(type, pvalue)), ]

results$Disease_or_Trait <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
library('ggplot2')
pdf("dSNP.eQTL.enrichment.pvalue.pdf", width=4, height=6)
results=subset(results, type!='gene')
p = ggplot(results, aes(x=Disease_or_Trait, y=-log10(pvalue), fill=type,ymax=max(-log10(pvalue))*1.2))
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity")
p = p + geom_hline(yintercept=-log10(0.05/1288), size=.5,linetype = 2)  ## Bonferroni correction, where 1288 is the number of disease/traits in GWAS
p = p + theme_bw() + ylab("-log10(p)") #+ scale_y_log10()
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8),legend.text=element_text(size=8), legend.justification='right', legend.position=c(1,1))
p = p + geom_text(aes(label=paste0(" (OR)")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5)
p = p + ggtitle(paste0(basename(getwd()), " -- dSNP enrichments"))
print(p)
dev.off()

quit('no')


#########################
### 5. post-normalization QC
#########################
Rscript $pipeline_path/modules/_normQC.R genes.fpkm.cufflinks.HCILB.uniq.xls QC.genes.fpkm.HCILB.uniq.pdf

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
