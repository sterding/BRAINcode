####################################
# Pipeline for CAGE data Analysis
# Authors: Xianjun DOng
# Email: xdong@rics.bwh.harvard.edu
# Date: 11/21/2014
# Version: 0.0
####################################
#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: $HOME/neurogen/pipeline/RNAseq/CAGE.pipeline.sh /data/neurogen/CAGE_PDBrainMap/rawfiles"
  exit
fi


########################
## 0. setting 
########################
pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

# project folders
input_dir=$1  # input_dir=/data/neurogen/CAGE_PDBrainMap/rawfiles

# create the subfolders (e.g. processed, for_display, results)
processed=$input_dir/../processed
[ -d $processed ] || mkdir $processed

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results 
[ -d $result_dir ] || mkdir $result_dir

########################
# 1. convert bam to bigwig
########################

cd $processed
for i in *bam; do
    bsub -J _bam2bigwig -oo _bam2bw.$i.log -eo _bam2bw.$i.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] -u $EMAIL -N _CAGE_bam2bigwig.sh $i;
done


