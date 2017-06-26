# ===============================================================
# convert input bed to narrowPeak with peak summit and signalValue annotated by $inputBigwig
# Author: Xianjun Dong
# Date: Jun 7, 2017
# ===============================================================

inputBed=$1  # at least a BED4
inputBigwig=$2

cat $inputBed | while read chr start end name rest
do
  l=$(expr $end - $start)
  bigWigSummary $inputBigwig $chr $start $end $l | awk -vchr=$chr -vstart=$start -vend=$end -vname=$name '{OFS="\t"; for(i=1;i<=NF;i++) if($i!="n/a" && $i>max) {imax=i;max=$i}}END{print chr, start, end, name, 0, ".", max, -1, -1, imax-1}'
done > $inputBed.narrowPeak
