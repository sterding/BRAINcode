###########################################
# script to generate coverage for combined samples
# Usage: $0 HC_SNDA 5 0.5
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
#!/bin/bash

group_lable=$1
N=$2  # at least N number of unique reads (if N is a int) or N rpm values (if N is float)
randN=$3 # randomly selected N samples

[ "$N" = "" ] && N=0.05;

n=`cat /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable | wc -l`
echo "total samples n = $n"

if [[ ! "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with reads cutoff"
    
    unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bw.gt5reads.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}reads.$group_lable.bed;
    sort -k4n,4n covered.${N}reads.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}reads.$group_lable.txt
fi

if [[ "$N" =~ "." ]]; then
    
    echo "["`date`"] calculating covered region with RPM cutoff"
    
    if [[ "$randN" = "" ]]; then
      unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > covered.${N}RPM.$group_lable.bed;
      sort -k4,4n covered.${N}RPM.$group_lable.bed | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > covered.${N}RPM.$group_lable.txt
    fi
    if [[ ! "$randN" = "" ]]; then
        TMPFILE=`mktemp random.covered.${N}RPM.$group_lable.XXXXXXXXXX.bed` || exit 1
        # reference: https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html
        # get_seeded_random()
        # {
        #   seed="$1"
        #   openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
        # }
    
        # To be fair, remove the samples with overall corerage (with RPM>0) less than 5% of the genome, e.g. those 18 SNDA samples at the low end of FigS1b.
        # People might be against this. TO MUTE THIS, just set lowcoveragesamples="NULL"
        lowcoveragesamples="HC_BN99-44_SNDA_5_rep1|HC_BN03-15_SNDA_5_rep1|ILB_BN06-25_SNDA_5_rep1|ILB_BN12-28_SNDA_5_rep1|ILB_BN99-50_SNDA_4_rep1|HC_BN02-24_SNDA_5_rep1|ILB_BN99-54_SNDA_4_rep1|HC_UK845_SNDA_6_rep1|HC_BN00-14_SNDA_4_rep1|HC_BN97-02_SNDA_4_rep1|HC_BN02-04_SNDA_5_rep1|HC_NZ-H83_SNDA_4_rep1|HC_BN04-52_SNDA_5_rep1|HC_BN03-41_SNDA_4_rep1|HC_BN08-64_SNDA_4_rep1|HC_BN07-11_SNDA_5_rep1|ILB_BN12-25_SNDA_5_rep1|HC_BN08-85_SNDA_5_rep1"
        unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -f /data/neurogen/rnaseq_PD/results/merged/samplelist.$group_lable | grep -vE $lowcoveragesamples | shuf -n $randN` | awk -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if(m>0) print $1,$2,$3,m}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > $TMPFILE
        sort -k4,4n $TMPFILE | awk '{OFS="\t"; print $0,$3-$2;}' | groupBy -g 4 -c 5 -o sum > $TMPFILE.txt
    fi  
fi



# for raw reads: at least N unique reads in at least M samples
#[[ ! "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"reads"M"samples:", s;}' > $group_lable.covered.${N}reads${M}samples.bed

# for normalized RPM: at least N RPM in at least M samples
#[[ "$N" =~ "." ]] && unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized.bedGraph | grep -E "$pattern"` | awk -vM=$M -vN=$N 'BEGIN{OFS="\t"; }{m=0; for(i=4;i<=NF;i++) if($i>=N) m++; if((M ~ /\./ && M<=1 && m/(NF-3)>=M) || (M !~ /\./ && M>=1 && m>=M)) {print $1,$2,$3}}' | mergeBed | awk -vM=$M -vN=$N 'BEGIN{OFS="\t";s=0;}{s+=($3-$2); print $0,1;}END{print "#"N"RPM"M"samples:", s;}' > $group_lable.covered.${N}RPM${M}samples.bed