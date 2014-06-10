#!/usr/bin/awk --posix -f
# awk script to parse pileup output from samtools to get depth of each variation in each sample
# Authos: Xianjun Dong
# Date: 2014--5-21
# Usage: awk -f _pileup2depth.awk input.pileup

## example of pileup ouput from samtools mpileup
#seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
#seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
#seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
#seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
#seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
#seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
#seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
#seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
#seq2 156 A 11  .$......+2AG.+2AG.+2AGGG    <975;:<<<<<
#seq3 200 A 20 ,,,,,..,.-4CACC.-4CACC....,.,,.^~. ==<<<<<<<<<<<::<;2<<
#seq2 200 A 235 .$.$.$.$.$.......................................................................................................................-1T......,...................................,........................................................^S.^S.^S.^S.^S.^S.^S.^S.^S.^S.^S.^S. FHDHD?HHGHGIFF<DBBFHBBGFHGFFFGG?FHFGFHH?G9GHGGHFFBBFHHGGFGHH9GFGGGG<FHHFGGFGHBHHB?@HFGDFG?:F7FGGGGDBG*FHGFFGFFGHF*GHGFHHGFHGFHEHFGJFHDGHGEGHGCFGGHGHC<HCFGFH?HHHHHGGJGIFGHGGF:HGHGHGFGCHHHHHFHHGFHHHHFHH:HHHGADFHHDDDADC@C@CCCC@C?@@@@@CCC<

BEGIN{
    OFS="\t";
}
{
    chr=$1; pos=$2; ref=toupper($3);
    printf "%s_%d_%s", chr, pos, ref;
    for(i=4;i<NF;i=i+3){
        rb=i+1; # read base field
        gsub("\\^.|[\\$\\-\\+]","",$rb);
        m=gsub("[,.]","",$rb);
        while(match($rb,"[1-9]+")>0) {
            n=substr($rb,RSTART, RLENGTH);
            gsub(n"[a-zA-Z]{"n"}","",$rb);
        }
        N["A"]=(ref=="A")?m:gsub("[Aa]","",$rb);
        N["G"]=(ref=="G")?m:gsub("[Gg]","",$rb);
        N["T"]=(ref=="T")?m:gsub("[Tt]","",$rb);
        N["C"]=(ref=="C")?m:gsub("[Cc]","",$rb);
        
        # all 4 nucleotides
        printf "\t%d:%d:%d:%d", N["A"], N["G"], N["T"], N["C"];     

    }
    printf "\n";
}
