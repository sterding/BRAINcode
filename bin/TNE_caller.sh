## General script to define TNE (rewritten from eRNA.define.sh)
# Author: Xianjun Dong
# Version: 2.0
# Date: Nov 15th, 2017
# Usage: TNE_caller.sh -h
# bsub -q long -n 2 -R 'rusage[mem=4000]' -J MSBB TNE_caller.sh -G MSBB 
# bsub -q long -n 2 -R 'rusage[mem=4000]' -J SKNSH TNE_caller.sh -G SKNSH
# bsub -q long -n 2 -R 'rusage[mem=4000]' -J TCPY2 TNE_caller.sh -G TCPY2
# bsub -q long -n 2 -R 'rusage[mem=4000]' -J BrainGVEX TNE_caller.sh -G BrainGVEX
# preqisition: 
#   TNE_caller.combine_bigwig.sh: generate trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph if it's not existed [optional]
#   TNE_caller.fit.Tx.noise.R: calculate RNAseq density with p=0.05
#   TNE_caller.consistency.R: calculate binominal pvalue from RPM in each sample

#!/bin/bash

################################################
# TNE definition:
# 1) density higher than the basal level,  
# 2) summit >0.05 RPM, --> p<0.05 comparing to the transcriptional noise
# 3) located in non-generic regions (e.g. 500bp away from any annotated exons),
# 4) at least 100bp in length,
# 5) don't contain any splicing sites (>10 reads in at least 5 samples) from Tophat
# 6) compute p-value for each TNE candidate given each sample's random background, then test the number of samples with p<0.05 with the binomial test, adjust p-value from binomial test with HB correction. Final TNEs have adjusted p-value <0.05.
################################################

Options=$@
Optnum=$#

## to generate the bigwig list file
## cd ~/eRNAseq
## mkdir MSBB; ls ~/neurogen/ROSMAP/MSBB/rnaseq_runoutput/*.bw | awk '{OFS="\t"; match($1,/runoutput\/(.*).Aligned/,a); print a[1], $1;}' > MSBB.bigwig.list 
## ls ~/neurogen/rnaseq_PD/run_output/CC_SK-N-SH-*/uniq/accepted_hits.normalized.bw | awk '{OFS="\t"; match($1,/run_output\/(.*)\/uniq/,a); print a[1], $1;}' > SKNSH.bigwig.list
## mkdir TCPY2; ls ~/neurogen/rnaseq_PD/run_output/*_TCPY_10_*/uniq/accepted_hits.normalized.bw | awk '{OFS="\t"; match($1,/run_output\/(.*)\/uniq/,a); print a[1], $1;}' > TCPY2/TCPY2.bigwig.list

function _usage()
{
  ###### U S A G E : Help and ERROR ######
  cat <<EOF
$0 $Options

   $*

      Usage: $0 <[options]>
      Options:
        -h   --help                 Show this message
        -G   --group_label=...      The label for samples group ($SAMPLE_GROUP)
        -B   --inputBG=...          The bedGraph file of merged samples ($inputBG)
        -g   --genome_size=...      genome chromosome sizes ($genome_size)
        -x   --excluded=...         BED file for exclusion region ($toExclude)
        -L   --list_bw_file=...     Tab-delimited file with two columns: sample ID and the path of bigwig file ($list_bw_file)
        -s   --splicing_site=...    BED file for splicing site +/-10bp ($splicing_site)
        -l   --length_min=...       minimal length of TNEs ($length_min)
EOF
}

[ $# = 0 ] && _usage "  >>>>>>>> no options given <<<<<<<< " && exit $?;

# ================================
# process the parameters
# ================================
# ref: https://stackoverflow.com/a/12523979

while getopts ':G:B:g:x:L:s:l:-h' OPTION ; do
  case "$OPTION" in
    h  ) _usage                                 ;;   
    G  ) SAMPLE_GROUP="$OPTARG"                 ;;
    B  ) inputBG="$OPTARG"                      ;;
    g  ) genome_size="$OPTARG"                  ;;
    x  ) excluded="$OPTARG"                     ;;
    L  ) list_bw_file="$OPTARG"                 ;;
    s  ) splicing_site="$OPTARG"                ;;
    l  ) length_min="$OPTARG"                   ;;
    -  ) optind=$OPTIND
         eval OPTION="\$$optind"
         OPTARG=$(echo $OPTION | cut -d'=' -f2)
         OPTION=$(echo $OPTION | cut -d'=' -f1)
         case $OPTION in
            --help  ) _usage                                           ;;   
            --group_label     ) SAMPLE_GROUP="$OPTARG"                 ;;
            --inputBG         ) inputBG="$OPTARG"                      ;;
            --genome_size     ) genome_size="$OPTARG"                  ;;
            --excluded        ) excluded="$OPTARG"                     ;;
            --list_bw_file    ) list_bw_file="$OPTARG"                 ;;
            --splicing_site   ) splicing_site="$OPTARG"                ;;
            --length_min      ) length_min="$OPTARG"                   ;;
             * )  _usage ">>>>>>>> invalid options <<<<<<<<" && exit $? ;;
         esac
       OPTIND=1
       shift
      ;;
    ? )  _usage ">>>>>>>> invalid options <<<<<<<<" && exit $? ;;
  esac
done

# SAMPLE_GROUP=MSBB
[ -z "$SAMPLE_GROUP" ] && _usage ">>>>>>>> invalid options <<<<<<<<" && exit $? ;
[ -z "$toExclude" ] && toExclude=/data/neurogen/external_download/externalData/toExclude.bed
[ -z "$length_min" ] && length_min=100
[ -z "$splicing_site" ] && splicing_site=/data/neurogen/rnaseq_PD/results/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed
[ -z "$genome_size" ] && genome_size=$GENOME/Annotation/Genes/ChromInfo.txt
[ -z "$inputBG" ] && inputBG=trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph
[ -z "$list_bw_file" ] && list_bw_file=$SAMPLE_GROUP.bigwig.list

[ -d $SAMPLE_GROUP ] || mkdir $SAMPLE_GROUP
cd $SAMPLE_GROUP

## quit if eRNA.bed already exists
[ -e eRNA.bed ] && echo "File eRNA.bed already exists. Quit!" && exit $?

[ -e $list_bw_file ] || ln -fs ../$list_bw_file $list_bw_file  # if not exist, then create a soft link in the current folder

echo "# step0: measure transcriptional noise in background genomic regions"
# =================

# make merged signal if not existed
[ -e $inputBG ] || bsub -q big -n 4 -J bams2combinedbg TNE_caller.combine_bigwig.v2.sh $list_bw_file $SAMPLE_GROUP && inputBG=trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph

# RNAseq signal distribution in the background region
[ -e transcriptional.noise.rpm.txt ] || bedtools random -seed 3 -g $genome_size -l 1 -n 1000000 | sortBed | intersectBed -a - -b $toExclude -v -sorted | intersectBed -a $inputBG -b - -sorted -u | cut -f4 > transcriptional.noise.rpm.txt
TNE_caller.fit.Tx.noise.R transcriptional.noise.rpm.txt
Dsig=`tail -n1 transcriptional.noise.rpm.pvalues.txt`  

echo "# step1: any regions with value > baseLevel (i.e. average sequencing depth)"
# =================
basalLevel=`tail -n1 $inputBG | cut -f2 -d'=' | cut -f1`
awk -vmin=$basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' $inputBG | mergeBed -c 5 -o max > eRNA.tmp1

echo "# step2: summit RPM >=Dsig (density with p<0.05)"
# =================
awk -vD=$Dsig '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' eRNA.tmp1 | mergeBed -d 100 -c 5 -o max > eRNA.tmp2

echo "# step3: located in non-generic regions (e.g. 500bp away from any annotated exons)"
# =================
intersectBed -a eRNA.tmp2 -b $toExclude -v > eRNA.tmp3

echo "# step4: length > $length_min"
# =================
awk -v len_min=$length_min '{OFS="\t"; if(($3-$2)>len_min) print $1,$2,$3,$1"_"$2"_"$3}' eRNA.tmp3 > eRNA.tmp4

echo "# step5: don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly)"
# =================
# cd ~/neurogen/rnaseq_PD/results/merged/denovo_assembly/
# cat cufflinks-cuffmerge/merged.bed trinity-cuffmerge/all_strand_spliced.chr.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' | bed12ToBed6 | awk '{OFS="\t"; print $1, $2-10,$2+10; print $1,$3-10,$3+10;}' | sortBed | uniq > trinitycufflinks.merged.splicingsites.flanking20nt.bed

# more than 10 splicing reads in at least 5 samples
# for i in  ~/neurogen/rnaseq_PD/run_output/*/junctions.bed; do awk '{OFS="\t"; if($5>10) { split($11,a,","); split($12,b,","); print $1,$2+a[1]-10,$2+a[1]+10; print $1,$2+b[2]-10,$2+b[2]+10}}' $i | sortBed | uniq; done | sort | uniq -c | awk '{OFS="\t"; if($1>5) print $2,$3,$4}' > ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed
 
intersectBed -a eRNA.tmp4 -b $splicing_site -v > eRNA.tmp5

echo "# step6: calculate the significance of eRNA"
# =================

#1: create random background regions (same number and same length distribution as TNEs) and calculate their signals
while read sample filename
do
    echo " - sample:" $sample;
    format=${filename##*.} # get extension
    
    if [[ $format =~ "bigwig|BIGWIG|bw|BW" ]];
    then
      bedtools shuffle -seed 123 -excl $toExclude -noOverlapping -i eRNA.tmp5 -g $genome_size | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bigWigAverageOverBed $filename stdin stdout | cut -f1,5 > $filename.$SAMPLE_GROUP.rdbg &
      bsub -q vshort -n 1 -J $sample "bigWigAverageOverBed $filename eRNA.tmp5 stdout | cut -f1,5 | sort -k1,1 -o $filename.$SAMPLE_GROUP.eRNA.meanRPM"
    elif [[ $format =~ "bam|BAM|cram|CRAM" ]]; 
    then
      bedtools shuffle -seed 123 -excl $toExclude -noOverlapping -i eRNA.tmp5 -g $genome_size | awk -vOFS="\t" '$4=$1"_"$2"_"$3;' | bedtools coverage -a stdin -b $filename -d | groupBy -g 4 -c 6 -o mean > $filename.$SAMPLE_GROUP.rdbg &
      bsub -q vshort -n 1 -J $sample "bedtools coverage -a eRNA.tmp5 -b $filename -d | groupBy -g 4 -c 6 -o mean | sort -k1,1 -o $filename.$SAMPLE_GROUP.eRNA.meanRPM"
    fi
    
done < $list_bw_file

#2. compute p-value for each TNE candidate in each sample's random background, then test the number of samples with p<0.05 with the binomial test, adjust p-value from binomial test with HB correction
TNE_caller.consistency.R $SAMPLE_GROUP $list_bw_file   # output are eRNA.tmp5.meanRPM.xls, eRNA.tmp5.pvalues.xls, and eRNA.tmp5.pvalues.adjusted.xls

#3. Select TNE with adjusted p <= 0.05: 
# bonferroni for the major groups
awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($4<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.bonferroni.bed
# FDR with method <=0.05
awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($5<=0.05) print a[1],a[2],a[3],$1}}' eRNA.tmp5.pvalues.adjusted.xls | sortBed > eRNA.fdr.bed

echo "# Final output files: eRNA.bed, eRNA.loci.txt, eRNA.meanRPM.xls"
# =================
# use bonferroni for major cell types and FDR for minor cell types.
ln -fs eRNA.fdr.bed eRNA.bed
# loci.txt required for fasteQTL
awk 'BEGIN{OFS="\t"; print "id","chr","s1","s2";}{print $4,$1,$2,$3;}' eRNA.bed > eRNA.loci.txt  
# meanRPM
paste eRNA.tmp5.pvalues.adjusted.xls eRNA.tmp5.meanRPM.xls | awk 'NR ==1 || $5<=0.05' | cut -f6- > eRNA.meanRPM.xls

echo "THE END!"





# ## expression level of the defined eRNAs on ALL (except the stranded) samples, not just the $SAMPLE_GROUP
# cut -f4 eRNA.bed | sort -k1,1 > eRNA.meanRPM.allSamples.xls
# TMPFILE=`mktemp /tmp/meanRPM.XXXXXXXXXX` || exit 1
# while read i
# do
#   echo $i;
#   bigWigAverageOverBed $i eRNA.bed stdout | sort -k1,1 | cut -f5 | paste eRNA.meanRPM.allSamples.xls - > $TMPFILE
#   cp $TMPFILE eRNA.meanRPM.allSamples.xls
# done < $list_bw_file_allSamples
# 
# echo "locus" `cat $list_bw_file_allSamples | sed 's/.*run_output\/\(.*\)\/uniq.*/\1/g' | rowsToCols stdin stdout` | sed 's/ /\t/g' | cat - $TMPFILE > eRNA.meanRPM.allSamples.xls