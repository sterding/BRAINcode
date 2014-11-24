###########################################
# script to combine bigwig, using trimmed mean of unionBedGraphs for samples from the same group
# Usage: $0 HC_SNDA
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
#!/bin/bash

group_lable=$1

[ "$group_lable" = "HC_TCPY" ]  && pattern="(HC|ND)_.*_TCPY_[2345]";
[ "$group_lable" = "HC_MCPY" ]  && pattern="(HC|ND)_.*_MCPY_[2345]";
[ "$group_lable" = "HC_SNDA" ]  && pattern="(HC|ND)_.*_SNDA_[2345]";
[ "$group_lable" = "ILB_SNDA" ] && pattern="ILB_.*_SNDA_[2345]";
[ "$group_lable" = "PD_SNDA" ]  && pattern="PD_.*_SNDA_[2345]";
[ "$group_lable" = "HCILB_SNDA" ]  && pattern="(HC|ND|ILB)_.*_SNDA_[2345]";

echo $group_lable, "$pattern";
ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern"
    
echo "["`date`"] computing trimmed mean of bedGraph"
#-------------------

unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern"` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); sum=0; for(i=a+1;i<=(c-a);i++) sum+=j[i];return sum/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > trimmedmean.uniq.normalized.$group_lable.bedGraph

# collapse continous regions with the same value by code: awk '{OFS="\t"; if(id!=$4 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}'

echo "["`date`"] bedgraph --> bigwig"
#-------------------

bedGraphToBigWig trimmedmean.uniq.normalized.$group_lable.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt trimmedmean.uniq.normalized.$group_lable.bw


## Note: other ways of averaging (e.g. mean, median etc.) of all samples [SAVE FOR REFERENCE]
## ================================================================

## version 1: mean of all samples

# multi-mapper
#unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_MCPY.mean.multi.normalized.bedGraph && bedGraphToBigWig HC_MCPY.mean.multi.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_MCPY.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > PD_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig PD_SNDA.mean.multi.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt PD_SNDA.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig HC_SNDA.mean.multi.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_SNDA.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > ILB_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.mean.multi.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt ILB_SNDA.mean.multi.normalized.bw

## uniq-mapper
#unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); if(s>0) print $1,$2,$3,s}' > HC_MCPY.mean.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.mean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_MCPY.mean.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); if(s>0) print $1,$2,$3,s}' > PD_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.mean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt PD_SNDA.mean.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); if(s>0) print $1,$2,$3,s}' > HC_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.mean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_SNDA.mean.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); if(s>0) print $1,$2,$3,s}' > ILB_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.mean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt ILB_SNDA.mean.uniq.normalized.bw
#
## version 2: median of all samples
#
#unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=median(S); if(tm>0) print $1,$2,$3,tm;}' > HC_MCPY.median.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.median.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_MCPY.median.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=median(S); if(tm>0) print $1,$2,$3,tm;}' > PD_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.median.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt PD_SNDA.median.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=median(S); if(tm>0) print $1,$2,$3,tm;}' > HC_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.median.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_SNDA.median.uniq.normalized.bw
#unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=median(S); if(tm>0) print $1,$2,$3,tm;}' > ILB_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.median.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt ILB_SNDA.median.uniq.normalized.bw


## version 3: trimmed mean (10%)
#
#unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "(HC|ND)_.*_MCPY_[2345]"` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' > HC_MCPY.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.trimmedmean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_MCPY.trimmedmean.uniq.normalized.bw
#
#unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "(HC|ND)_.*_TCPY_[2345]"` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' > HC_TCPY.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig HC_TCPY.trimmedmean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_TCPY.trimmedmean.uniq.normalized.bw
#
#unionBedGraphs -i `ls ../../run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "(HC|ND)_.*_SNDA_[2345]"` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' > HC_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.trimmedmean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt HC_SNDA.trimmedmean.uniq.normalized.bw
#
#unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/accepted_hits.normalized2.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' > PD_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.trimmedmean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt PD_SNDA.trimmedmean.uniq.normalized.bw
#
#unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/accepted_hits.normalized2.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm>0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' > ILB_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.trimmedmean.uniq.normalized.bedGraph $GENOME/Annotation/Genes/ChromInfo.txt ILB_SNDA.trimmedmean.uniq.normalized.bw


