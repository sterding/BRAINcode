#!/bin/awk -f
# awk script to filter noise by sliding a window and taking mean per window
# Authos: Xianjun Dong
# Date: 2015-06-23
# Usage: _filter.awk values.txt
# bigWigSummary externalData/RNAseq/RNAseq.BRAINCODE.HCILB_SNDA.bigwig chr10 101130293 101131543 104 | _filter.awk -vW=5 -vtype=median
BEGIN{
  if(W=="") W=5; 
  if(type=="") type="median";
  half=(W-1)/2;
} 
{
  for(i=1;i<=NF;i++) {
    array[half+1]=$i
    for(k=1;k<=half;k++){
      array[half+1-k]=(i<=k)?$1:$(i-k);
      array[half+1+k]=((i+k)>NF)?$NF:$(i+k);
    }
    if(type=="median") {asort(array); x=array[half+1];}
    if(type=="mean") {x=0; for(j=1;j<=W;j++) x+=array[j]; x=x/W;}
    printf("%s\t", x);
  }
  print "";
}