#!/usr/bin/awk -f
# awk script to calculate normalized CpG score for input DNA sequences
# Authos: Xianjun Dong
# Date: 2015-03-30
# Usage: getNormalizedCpGscore.awk input.fa.tab
# the input.fa.tab is a text file with two columns: ID and sequence. It can be extracted by tools like "bedtools getfasta -name -tab -fi "

# The calculation is similar as https://github.com/sterding/code_umass/blob/master/piRNA2/get_normalized.CpG.content.pl
# Normalized CpG fraction was computed as (observed CpG)/(expected CpG), where expected CpG was calculated as (GC content/2)^2.

{
    OFS="\t";
    id = $1;
    s = tolower($2);
    
    
    n1 = gsub("gc","gc",s);
    n2 = length($s);
    n3 = gsub("g","g",s) + gsub("c","c",s);

    print id, 4*n2*n1/(n3**2);
}