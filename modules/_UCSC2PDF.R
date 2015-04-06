# Script: converting UCSC region (in bed format) to PDF files
# Usage: Rscript _UCSC2PDF.R bedfile
source('~/neurogen/pipeline/RNAseq/bin/lib.R')

args <- commandArgs(TRUE)
bedfile=args[1]

zoomoutX=2

# example of controling individual track
#theURL="http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&wgRna=hide&cpgIslandExt=pack&ensGene=hide&mrna=hide&intronEst=hide&mgcGenes=hide&hgt.psOutput=on&cons44way=hide&snp130=hide&snpArray=hide&wgEncodeReg=hide&pix=1000&refGene=pack&knownGene=hide&rmsk=hide"
# example of controling via session
theURL="http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sterding&hgS_otherUserSessionName=hg19_PD&hgt.psOutput=on&pix=1000"
#library(multicore)

# read regions
toPlot=read.table(bedfile, header=F)

for(i in 1:nrow(toPlot)){
    length=abs(toPlot[i,3]-toPlot[i,2])
    chr=as.character(toPlot[i,1])
    start=toPlot[i,2]-round(length*(zoomoutX-1)/2)
    end=toPlot[i,3]+round(length*(zoomoutX-1)/2)
    filename=paste0("/tmp/eRNAucsc.", chr,"_",toPlot[i,2],"_",toPlot[i,3],".pdf")
    
    if(file.exists(filename)) next;
        
    screenshotUCSC(theURL, "", chr, start, end, filename)

    ## convert to PNG
    #command=paste("convert -density 120", filename, gsub(".pdf",".png", filename))
    #cat(command,"\n")
    #try(system(command))
    #
    ### scp to zlab
    #command=paste("chmod 644", gsub(".pdf",".p*", filename), "; scp", gsub(".pdf",".p*", filename), "panda.dipr.partners.org:~/public_html/rnaseq_PD/version2/eRNA/ucsc")
    #cat(command,"\n")
    #try(system(command))

    # anti-robot version
    #Sys.sleep(5) # Program-driven use of this software is limited to a maximum of one hit every 15 seconds and no more than 5,000 hits per day.
    print(paste(i, filename));
}