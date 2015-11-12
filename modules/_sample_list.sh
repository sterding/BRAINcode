###########################################
# script to generate sample list per group
# Usage: $0 HC_SNDA
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
#!/bin/bash

group_lable=$1

# dopamine neurons (HC+ILB, SNDA)
[ "$group_lable" = "HCILB_SNDA" ]  && pattern="(HC|ILB)_.*_SNDA_[^1]_[^s]*/";  # excep batch1 and stranded libs
# pyramidal neurons (HC, TCPY+MCPY)
[ "$group_lable" = "HC_PY" ]  && pattern="HC_.*_[TM]CPY_[^1]";
# non-neuronal (HC, PBMC+FB)
[ "$group_lable" = "HC_nonNeuron" ]  && pattern="HC_.*_(FB|PBMC)_[^u]*/";
# neuronal (HC+ILB, SNDA+TCPY+MCPY)
[ "$group_lable" = "HC_Neuron" ]  && pattern="(HC|ILB)_.*_(SNDA|TCPY|MCPY)_[^1]_[^s]*/";

[ "$group_lable" = "HC_TCPY" ]  && pattern="HC_.*_TCPY_[^1]";
[ "$group_lable" = "HC_MCPY" ]  && pattern="HC_.*_MCPY_[^1]";
[ "$group_lable" = "HC_SNDA" ]  && pattern="HC_.*_SNDA_[^1]_[^s]*/";
[ "$group_lable" = "ILB_SNDA" ] && pattern="ILB_.*_SNDA_[^1]";
[ "$group_lable" = "PD_SNDA" ]  && pattern="PD_.*_SNDA_[^1]";
# control
[ "$group_lable" = "HC_SN" ]  && pattern="HC_.*_SN_[^u]*/";
[ "$group_lable" = "HC_SNDAstranded" ]  && pattern="HC_.*_SNDA_.*stranded";
[ "$group_lable" = "HC_PBMC" ]  && pattern="HC_.*_PBMC_[^u]*/";
[ "$group_lable" = "HC_FB" ]  && pattern="HC_.*_FB_";

ls ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized2.bedGraph | grep -E "$pattern" | sed 's/.*run_output\/\(.*\)\/uniq.*/\1/g' 