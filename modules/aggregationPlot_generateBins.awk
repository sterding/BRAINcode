# ===============================================================
# awk script to convert a transcript region to 81 bins (40bins around TSS, 40bins around TTS, and 1 bin for the rest of gene body)
# Author: Xianjun Dong
# Date: Mon Feb 10
# the input file is a gtf/bed file, each row is a transcript with columns chr/TSS/TTS/strand/ID
# the output file is a bed file for all bins, and it will be used for input of bigWigAverageOverBed (fron Jim kent)
# Usage: cat input.bed12 | awk -v option=mRNA -f bigWigAverageOverBed_generate_81bins.awk
# ===============================================================
#!/bin/awk -f

# fixed-size bins for flanking region, proportional-sized bins for gene body region (without intron); input is bed12 format
{
        chr=$1; start=$2; end=$3; name=$4; strand=$6; exonCount=$10; exonSizes=$11; exonStarts=$12;
        if((end-start)<41) next;

        #print option;

        if(option=="mRNA" || option=="cDNA" || option=="RNA" || option=="transcriptome")
        {
            FLANKING=1000;
            if(start<FLANKING) FLANKING=start;
            #print FLANKING;
            if(FLANKING<20) next;

            BINSIZE=FLANKING/20;  # bin size for flanking region

            # sum of exons length
            exonlengthSum=0;
            split(exonStarts, b, ","); split(exonSizes, c, ",");
            for(i=1;i<=length(c);i++) exonlengthSum=exonlengthSum+c[i];
            
            BINsize=exonlengthSum/41; # bin size for exonic region
        

            if(strand=="+"){
                # upstream (20 bins)
                for(i=0;i<20;i++) print chr, start-FLANKING+int(BINSIZE*i+0.5),$2-FLANKING+int(BINSIZE*(i+1)+0.5), name".bin"(i+1);

                # gene body (41 bins)
                split(exonStarts, b, ","); split(exonSizes, c, ","); j=0;
                for(i=1;i<exonCount;i++) {
                    binstart=start+b[i]; exonEnd=start+b[i]+c[i];
                    #print exonlengthSum, BINsize,binstart, exonEnd, end, j, i, exonCount;
                    while((binstart+int(BINsize+0.5)) <= exonEnd)
                    {
                        if(j>=40) break;
                        print chr, binstart, binstart+int(BINsize+0.5), name".bin"(j+21);
                        binstart=binstart+int(BINsize+0.5);
                        j++;
                    }
                }
                #print "split the most 3\' exon into the rest of (41-j) bins", j, i, exonCount, b[exonCount], length(b)
                binsize=c[exonCount]/(41-j)
                for(k=j;k<41;k++) print chr, start+b[exonCount]+int(binsize*(k-j)+0.5),start+b[exonCount]+int(binsize*(k-j+1)+0.5), name".bin"(k+21);

                #downstream
                #print "split the downstream into the 20 bins"
                for(i=0;i<20;i++) print chr, end+int(BINSIZE*i+0.5), end+int(BINSIZE*(i+1)+0.5), name".bin"(i+62);
            }
            if(strand=="-"){
                # upstream (20 bins)
                for(i=0;i<20;i++) print chr, end+FLANKING-int(BINSIZE*(i+1)+0.5),end+FLANKING-int(BINSIZE*i+0.5), name".bin"(i+1);

                # gene body (41 bins)
                j=0;
                for(i=exonCount;i>1;i--) {
                    exonStart=start+b[i]; binend=start+b[i]+c[i];
                    while((binend-BINsize) >= exonStart)
                    {
                        if(j>=40) break;
                        print chr, int(binend-BINsize), int(binend), name".bin"(j+21);
                        binend=binend-BINsize;
                        j++;
                    }
                }
                #print "split the most 3\' exon into the rest of (41-j) bins", j, i, exonCount, b[exonCount], length(b)
                binsize=c[1]/(41-j)
                for(k=j;k<41;k++) print chr, start+b[1]+c[1]-int(binsize*(k-j+1)+0.5),start+b[1]+c[1]-int(binsize*(k-j)+0.5), name".bin"(k+21);

                #print "downstream (20 bins)"
                for(i=0;i<20;i++) print chr, start-int(BINSIZE*(i+1)+0.5), start-int(BINSIZE*i+0.5), name".bin"(i+62);
            }
        }
        
        if(option=="genome" || option=="" || option=="DNA")
        {
            # flanking=min(1000, start, length/2)
            FLANKING=1000;
            if(start<FLANKING) FLANKING=start;
            if((end-start-1)<2*FLANKING) FLANKING=int((end-start-1)/2);
            #print FLANKING;

            if(FLANKING<20) next;
            BINSIZE=FLANKING/20;

            if(strand=="+"){
                for(i=0;i<40;i++) print chr, start-FLANKING+int(BINSIZE*i),start-FLANKING+int(BINSIZE*(i+1)), name".bin"(i+1);
                print chr, start-FLANKING+int(BINSIZE*40), end-FLANKING, name".bin41";
                for(i=0;i<40;i++) print chr, end-FLANKING+int(BINSIZE*i), end-FLANKING+int(BINSIZE*(i+1)), name".bin"(i+42);
            }
            if(strand=="-"){
                for(i=0;i<40;i++) print chr, end+FLANKING-int(BINSIZE*(i+1)), end+FLANKING-int(BINSIZE*i), name".bin"(i+1);
                print chr, start+FLANKING, end+FLANKING-int(BINSIZE*40), name".bin41";
                for(i=0;i<40;i++) print chr, start+FLANKING-int(BINSIZE*(i+1)), start+FLANKING-int(BINSIZE*i), name".bin"(i+42);
            }
        }
}

# Usage: cut -f1-4,6 /home/dongx/nearline/genomes/mm9/Annotation/Genes/NCBIM37.biomart67.transcripts.cpg.tab | awk -v option=genome -f bigWigAverageOverBed_generate_81bins.awk > ../data/NCBIM37.biomart67.transcripts.cpg.tab.81bins &
## fixed-size bins for flanking region, proportional-sized bins for gene body region (81bins)
#{
#        # flanking=min(1000, start, length/2)
#        FLANKING=1000;
#        if($2<FLANKING) FLANKING=$2;
#        if(($3-$2-1)<2*FLANKING) FLANKING=int(($3-$2-1)/2);
#        #print FLANKING;
#
#        if(FLANKING<20) next;
#
#        BINSIZE=FLANKING/20;
#
#        if($4=="+"){
#            for(i=0;i<40;i++) print $1, $2-FLANKING+int(BINSIZE*i),$2-FLANKING+int(BINSIZE*(i+1)), $5".bin"(i+1);
#            print $1, $2-FLANKING+int(BINSIZE*40), $3-FLANKING, $5.".bin41";
#            for(i=0;i<40;i++) print $1, $3-FLANKING+int(BINSIZE*i), $3-FLANKING+int(BINSIZE*(i+1)), $5.".bin"(i+42);
#        }
#        if($4=="-"){
#            for(i=0;i<40;i++) print $1, $3+FLANKING-int(BINSIZE*(i+1)),$3+FLANKING-int(BINSIZE*i), $5".bin"(i+1);
#            print $1, $2+FLANKING, $3+FLANKING-int(BINSIZE*40), $5.".bin41";
#            for(i=0;i<40;i++) print $1, $2+FLANKING-int(BINSIZE*(i+1)), $2+FLANKING-int(BINSIZE*i), $5.".bin"(i+42);
#        }
#}

# fixed-size bins for flanking region, proportional-sized bins for gene body region (81bins)
#{
#        FLANKING=2000;
#        if($4=="+"){
#            if($2<FLANKING) FLANKING=$2; BINSIZE=int(FLANKING/20);
#            for(i=0;i<20;i++) print $1, $2-FLANKING+BINSIZE*i,$2-FLANKING+BINSIZE*(i+1), $5".bin"(i+1);
#            BINSIZE=int(($3-$2)/41);
#            for(i=0;i<40;i++) print $1, $2+BINSIZE*i,$2+BINSIZE*(i+1), $5".bin"(i+21);
#            print $1, $2+BINSIZE*40,$3, $5".bin61";
#            BINSIZE=int(FLANKING/20);
#            for(i=0;i<20;i++) print $1, $3+BINSIZE*i,$3+BINSIZE*(i+1), $5".bin"(i+62);
#        }
#        if($4=="-"){
#            BINSIZE=int(FLANKING/20);
#            for(i=0;i<20;i++) print $1, $3+FLANKING-BINSIZE*(i+1),$3+FLANKING-BINSIZE*i, $5".bin"(i+1);
#            BINSIZE=int(($3-$2)/41);
#            for(i=0;i<40;i++) print $1, $3-BINSIZE*(i+1),$3-BINSIZE*i, $5".bin"(i+21);
#            print $1, $2,$3-BINSIZE*40, $5".bin61";
#            if($2<FLANKING) FLANKING=$2; BINSIZE=int(FLANKING/20);
#            for(i=0;i<20;i++) print $1, $2-BINSIZE*(i+1), $2-BINSIZE*i, $5".bin"(i+62);
#        }
#}

# fixed-size bins (prons: each bin has same size; cons: for smaller genes, this can induce bias)
#{
#        FLANKING=2000;
#        BINSIZE=int(2*FLANKING/40);
#        if($4=="+"){
#            if($2<FLANKING) FLANKING=$2;
#            for(i=0;i<40;i++) print $1, $2-FLANKING+BINSIZE*i,$2-FLANKING+BINSIZE*(i+1)-1, $5".bin"(i+1);
#            print $1, $2+FLANKING, $3-FLANKING-1, $5.".bin41";
#            for(i=0;i<40;i++) print $1, $3-FLANKING+BINSIZE*i, $3-FLANKING+BINSIZE*(i+1)-1, $5.".bin"(i+42);
#        }
#        if($4=="-"){
#            if($3<FLANKING) FLANKING=$3;
#            for(i=0;i<40;i++) print $1, $3+FLANKING-BINSIZE*(i+1)+1,$3+FLANKING-BINSIZE*i, $5".bin"(i+1);
#            print $1, $2+FLANKING+1, $3-FLANKING, $5.".bin41";
#            for(i=0;i<40;i++) print $1, $2+FLANKING-BINSIZE*(i+1)+1, $2+FLANKING-BINSIZE*i, $5.".bin"(i+42);
#        }
#}