#!/usr/bin/env Rscript
###########################################
# script to generate sample list per group
# Usage: $0 HC_SNDA
# for i in HCILB_SNDA HC_PY HC_nonNeuron HC_Neuron HC_MCPY HC_TCPY HC_SNDA ILB_SNDA HC_PBMC HC_FB PD_SNDA HC_SNDAstranded; do echo $i; _sample_list.R $i samplelist.$i; done
# Author: Xianjun Dong
# Version: 1.0
# Date: 2014-Oct-22
###########################################
args<-commandArgs(TRUE)
group_lable=args[1]
output=args[2]

if(length(args)!=2) stop("Usage: _sample_list.R HCILB_SNDA output.file")

pattern = switch(group_lable, 
                 # major groups
                 HCILB_SNDA = "(HC|ILB)_.*_SNDA", # dopamine neurons (HC+ILB, SNDA)
                 HC_PY = "HC_.*_[TM]CPY", # pyramidal neurons (HC, TCPY+MCPY)
                 HC_nonNeuron = "HC_.*_(FB|PBMC)_", # non-neuronal (HC, PBMC+FB)
                 # minor groups
                 HC_Neuron = "(HC|ILB)_.*_(SNDA|TCPY|MCPY)", # neuronal (HC+ILB, SNDA+TCPY+MCPY)
                 HC_TCPY = "HC_.*_TCPY",
                 HC_MCPY = "HC_.*_MCPY",
                 HC_SNDA = "HC_.*_SNDA_",
                 ILB_SNDA = "ILB_.*_SNDA",
                 HC_PBMC = "HC_.*_PBMC_",
                 HC_FB = "HC_.*_FB_",
                 PD_SNDA = "PD_.*_SNDA",
                 # controls
                 HC_SN = "HC_.*_SN_",
                 HC_SNDAstranded = "HC_.*_SNDA_.*stranded",
                 stop("Unrecognized group label!")
                );

# covaraice table URL
covURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"
if (!require("RCurl",quietly =T)) {
  install.packages("RCurl", quiet =T, dependencies = TRUE)
  library(RCurl,quietly = T)
}

cov=read.delim(textConnection(getURL(covURL)), stringsAsFactors =F)

cov = subset(cov, grepl(pattern, sampleName))
if(group_lable!='HC_SNDAstranded' & group_lable!='PD_SNDA') cov=subset(cov, BRAINCODE.final.selection==1)

write.table(cov$sampleName, file = output, quote=F, col.names = F, row.names = F)