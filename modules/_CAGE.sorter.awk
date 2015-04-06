#!//bin/awk -f
# awk script to sort CAGE data according to barcode, remove the following EcoP15I site CAGCAG, the first 5'nt G (if any), and 3' linker (if any). Save the remained reads into individual files (named with barcode), if it's >=10 nt in length.
# hopefully to minic the tagDust: tagdust -1 B:AGA,CTT,GAT -2 S:CAGCAG -3 G:G -4 R:N -5 P:TCGTATGCCGTCTTCTGCTTG $i -o $i -ref $GENOME/Annotation/Genes/U13369.1.fa -fe 2
# Authos: Xianjun Dong
# Date: 2015-1-28
# Usage: _CAGE.sorter.awk -v BARCODE="AGA|CTT|GAT" input.fastq
#

{
    OFS="\n";
    s[NR%4]=$0;
    if(NR%4==0) {
        if(s[1]!~/:Y:/ && match(s[2], "^("BARCODE")CAGCAGG?")) {
            bc=substr(s[2],1,3);
            s[2]=substr(s[2],RLENGTH+1,length(s[2]));
            s[0]=substr(s[0],RLENGTH+1,length(s[0]));
            if(match(s[2],"TCGTATG.*") || match(s[2],"TCGTAT$") || match(s[2],"TCGTA$") || match(s[2],"TCGT$") || match(s[2],"TCG$") || match(s[2],"TC$")) {
                s[2]=substr(s[2],1,length(s[2])-RLENGTH);
                s[0]=substr(s[0],1,length(s[0])-RLENGTH);
            }
            if(length(s[2])>=10) print s[1], s[2], s[3], s[0] >> "sorted."bc".fq";
        }
        else print s[1], s[2], s[3], s[0];
    }
}

## equailent perl code
#cat test.fq | perl -ne '$h=$_; $s=<>; $s2=<>; $s3=<>;if($h!~/:Y:/ && $s=~s/^((AGA|CTT|AAT)CAGCAGG?)//) {$bc=$2; $s3=substr($s3, length($1), length($s3)); $s3=sprintf("%s\n", substr($s3,0, length($s3)-length($1)-1)) if($s=~s/(TCGTATG.*)//); print $h,$s,$s2,$s3;}'

