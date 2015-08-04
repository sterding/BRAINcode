#!/bin/awk -f
# awk script to detect the bimodel (peak-valley-peak) intervals from a series of values
# the output is in narrowPeak format, where col5=length_of_interval, col6=., col7=height-of-left-peak, col8=height-of-valley, col9=height-of-right-peak, col10=length-from-valley-to-left-peak
# Authos: Xianjun Dong
# Date: 2015-06-23
# Usage: _wave_detector.awk -vchr=chr10 -vstart=101130293 -vend=101131543 values.txt
# bigWigSummary externalData/RNAseq/RNAseq.BRAINCODE.HCILB_SNDA.bigwig chr10 101130293 101131543 104 | awk '{for(i=3;i<=NF-2;i++) {print ($(i-2)+$(i-1)+$i+$(i+1)+$(i+2))/5}}' |tr '\n' '\t'| _wave_detector.awk -vchr=chr10 -vstart=101130293 -vend=101131543 


BEGIN{
  max1=0;min=999999;max2=0;
  OFS="\t";
}
{
  ind=(end-start)/NF;
  for(i=1;i<=NF;i++) {
    #print $i, "---", imax1, max1, imin, min, imax2, max2, "---";
    if($i>max1 && min==999999) {max1=$i; imax1=i;}
    else if(max1>0 && max2==0 && $i<min) {min=$i; imin=i;}
    else if(max1>0 && min!=999999 && $i>=max2) {max2=$i; imax2=i;}
    else if(max2>0 && $i<max2) {
      print chr, start+int((imax1-0.5)*ind), start+int((imax2-0.5)*ind), ".", int(ind*(imax2-imax1)), ".", max1, min, max2, int((imin-imax1)*ind);
      max1=max2; imax1=imax2; min=$i; imin=i; max2=0;
    }
  }
}
END{
  if(max2!=0) print chr, start+int((imax1-0.5)*ind), start+int((imax2-0.5)*ind), ".", int(ind*(imax2-imax1)), ".", max1, min, max2, int((imin-imax1)*ind);
}


 
 
 