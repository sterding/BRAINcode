#!/usr/bin/awk -f
# awk script to detect variation (SNP, indel etc.) in a input SAM file
# Authos: Xianjun Dong
# Date: 2013-11-11
# Usage: awk -f _sam2variation.awk input.sam

{
        OFS="\t";
        if($0~/^@/) next;
        
        readlength=length($10);
        
        # detecting SNP
        if($6==readlength"M"){
            if($0~/[ \t]MD:Z:{readlength}[ \t\n]/) next;
            if($0!~/[ \t]MD:Z:[0-9]+([A-Z][0-9]+)*[ \t\n]/) next;  # only take case with variation, not indel
            for(i=6;i<=NF;i++)
                if($i~/^MD:Z:/) {
                    MDline=substr($i, 6);
                    where=match(MDline, "^[0-9]+");
                    L=substr(MDline, 1, RLENGTH);
                    MD=substr(MDline, RLENGTH+1);
                    while(match(MD, "^[A-Z][0-9]+")) {
                        v=substr(MD, 1,1)
                        # get the positon of variation
                        chr=$3;start=$4+L;
                        # get the sequence in reads
                        s=substr($10, L+1, 1)
                        # print out the variation
                        #print $1, $10, $i, L, chr":"start, v":"s;
                        if(s!="N") print chr":"start-1"-"start, v":"s, $1, L+1;
                        L=L+substr(MD, 2, RLENGTH-1)+1;
                        MD=substr(MD, RLENGTH+1);
                    }
                    next;
                }
        }
        
        # detecting insertion
        if($6~/I/){
               #print;
               #TOADD:
        }
        # detecting deletion
        if($6~/D/){
               #print;
               #TOADD:
        }
}