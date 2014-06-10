###########################################
# R script for running pathway analysis using SPIA and DAVID
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 1/27/2014
# version: 1.0
# Note: 
###########################################

## TODO: use the updated source code of SPIA (from Junko) to avoid the error message


calcPathway <- function(file, bg, out="out", sp="hsa")
{

   library(SPIA)

   # store data to vectors
   data    <- read.table(file, sep="\t")

   allGIs <- as.character(read.table(bg)[,1])
   degGIs <- as.vector(data[,4])
   names(degGIs) <- as.vector(data[,1])
   res <- spia(de=degGIs, all=allGIs, organism=sp, nB=2000, plots=FALSE, beta=NULL, combine="fisher")

   # make text file and plot
   write.table(res, paste(out, ".txt", sep=""), sep="\t", append=F, quote=F, row.names=F, col.names=F)

   return(res)
}

# GO/Pathway analysis by DAVID
# http://bioconductor.org/packages/release/bioc/html/DAVIDQuery.html
