#######################################
#  script to calculate reads count in $1 input per interval in $2
# Author: Xianjun Dong
# Date: 2013-12-14
# Version: 0.0
#######################################

inputbam=$1  # bam format
intervals=$2  # bed format
outputfile=$3
    
coverageBed -abam -split -counts -a $inputbam -b $intervals | sort -k4,4 > $outputfile



