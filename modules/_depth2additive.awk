#!/usr/bin/gawk --posix -f
# awk script to compute additive genotype mode by parsing the depth output from awk script _pileup2depth.awk
# Authos: Xianjun Dong
# Date: 2014--5-21
# Usage: awk -f _depth2additive.awk input.pileupoutput

## input format:
#chrM_35_G	0:273:0:0	0:265:0:0	0:242:0:0	0:278:0:0
#chrM_54_G	0:291:0:0	0:284:0:0	0:241:0:0	0:294:0:0
#chrM_55_T	0:0:290:2	0:0:285:0	0:1:242:0	1:0:294:0
#chrM_56_A	293:0:0:0	286:0:0:0	245:0:0:0	296:0:0:0
#chrM_60_T	0:0:297:0	0:0:290:0	0:0:239:0	2:0:298:0

BEGIN{
    OFS="\t";
    I["A"]=1; I["G"]=2; I["T"]=3; I["C"]=4;
    minDP=50;
    minSF=0.5;
    minAF=0.1; # alternative allele frequency: at least 10% of total reads
}
{
    split($1, ID, "_");
    ref=ID[3];
    n=0; # number of samples with at least $minDP reads sequenced at the position
    for(i=2;i<=NF;i++){
        split($i, N, ":");
        sum=N[1]+N[2]+N[3]+N[4];
        if(sum >= minDP) n++;
        asort(N, sN);
        #asorti(N, sNi);
        $i=(sum==0)?"NA":((sN[3] > sum*minAF && sN[3]>5 )?1:((N[I[ref]]==sN[4])?0:2));
        #$i=(sum<minDP)?"NA":N[I[ref]]/sum;
    }
    if(n >= minSF*(NF-1)) print;
}