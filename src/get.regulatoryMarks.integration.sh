#============================================================
# Script to download bigwig files for various marks
# Usage: bash $pipeline_path/src/get.regulatoryMarks.integration.sh H3K4me1
# bsub -q big -n 5 -M 40000  bash $pipeline_path/src/get.regulatoryMarks.integration.sh H3K4me1
# bsub -q big -n 5 -M 40000  bash $pipeline_path/src/get.regulatoryMarks.integration.sh H3K27ac
#============================================================
## 
#!/bin/bash

cd ~/neurogen/external_download/roadmap

mark=$1

# echo "# download the data from Roadmap"
# echo "##############################"
# # ref: https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15
# #for i in E071 E074 E068 E069 E072 E067 E073 E070 E082 E081;  # brain only
# > url.$mark
# for i in E017 E002 E008 E001 E015 E014 E016 E003 E024 E020 E019 E018 E021 E022 E007 E009 E010 E013 E012 E011 E004 E005 E006 E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046 E030 E026 E049 E025 E023 E052 E055 E056 E059 E061 E057 E058 E028 E027 E054 E053 E112 E093 E071 E074 E068 E069 E072 E067 E073 E070 E082 E081 E063 E100 E108 E107 E089 E090 E083 E104 E095 E105 E065 E078 E076 E103 E111 E092 E085 E084 E109 E106 E075 E101 E102 E110 E077 E079 E094 E099 E086 E088 E097 E087 E080 E091 E066 E098 E096 E113 E114 E115 E116 E117 E118 E119 E120 E121 E122 E123 E124 E125 E126 E127 E128 E129;
# do 
#   echo http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/$i-$mark.pval.signal.bigwig >> url.$mark
# done
# parallel -a url.$mark axel -n 5 

# number of samples for the mark
N=`ls *$mark.pval.signal.bigwig | wc -l`
echo "# In total, $N datasets downloaded from Roadmap for mark $mark"

echo "# merge the signal into one track"
echo "##############################"
bigWigMerge *$mark.pval.signal.bigwig merged.$mark.pval.signal.bg
bedGraphToBigWig merged.$mark.pval.signal.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size merged.$mark.pval.signal.bigwig

echo "# call peaks from merged signal"
echo "##############################"
module load macs2/2.1
macs2 bdgpeakcall -i merged.$mark.pval.signal.bg --min-length 100 -o merged.$mark.pval.signal.p0.01.halfsamples.peaks --cutoff $N  # at least half samples with p=0.01 (2*$N/2)

echo "# job submitted!!"