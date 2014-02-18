plotClassification <- function(all, subtype)
{
    # all=genes;  subtype="protein_coding"
    subset=all;
    if(subtype!='all') subset=all[grep(subtype, rownames(all)),]

    subset=subset[!is.na(rowSums(subset)),] # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
    subset=subset[rowSums(subset)>0,] # Optional: remove genes with all zero (or unexpressed)
    subset=subset*1e6 # TEMPERAL USE: "uniq" dataset was not properly normalized

    # only the gene body
    gb=subset[,21:61]
    gb=gb[rowSums(gb)>0,] # Optional: remove genes with all zero (or unexpressed)
    
    # only the protein-coding transcripts (from protein-coding genes) [test]
    gb=gb[grep("\\.protein_coding",rownames(gb)),]
    
    d <- dist(gb, method = "euclidean")
    hc <- hclust(d, method="average")
    #plot(hc) # display dendogram
    image(d[[hc$order]])
}

getAGGREGATION <- function(all, subtype,  TRIM=0, meanmedian='mean')
{
    subset=all;
    if(subtype!='all') subset=all[grep(subtype, rownames(all)),]

    subset=subset[!is.na(rowSums(subset)),] # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
    subset=subset[rowSums(subset)>0,] # Optional: remove genes with all zero (or unexpressed)
    subset=subset*1e6 # TEMPERAL USE: "uniq" dataset was not properly normalized

    # all genes
    p1=subset; n1=nrow(p1)
    # remove gene body bin, and 1% outlier
    if(meanmedian=='mean') p1=apply(p1, 2, function(x) mean(x, trim=TRIM));
    if(meanmedian=='median') p1=apply(p1, 2, median);
    #p1[41]=NA
    
    return(list(p1, nrow(subset)))
}

# aggregation plot for non-strand data
draw.plot <- function(p,ylim=range(p,na.rm=T), main="", xax=T, legend="", ylab="Mean signal", yax=T)
{
             #plot
             par(mar=c(2,2,1,1), lwd=2/3)
             plot(p, ylim=ylim, ylab='', xlab='', xaxt='n',yaxt='n', main='', type='n')
             if(xax==T) axis(1, at=c(1,20,61,81), labels=c("-1k bp","5'","3'","+1k bp"), tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1.5)
             if(yax==T) axis(4, las=1, tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1)
             mtext(ylab, side=2, line=1, cex=1)
             mtext(main, side=3, line=1, cex=1, padj=0)
             lines(c(1,81),c(0,0))
             #lines(c(42,81),c(0,0))
             abline(v=c(20,61), lty=2, col='darkgray')
             polygon(c(1,1:81,81),c(0,p[1:81],0), col='red', border=NA)
             if(legend!="") legend("topleft",legend, bty='n', cex=1.5)
}


# Here is an R script wrote by Aaron Statham which saves UCSC to pdfs -
# you can choose which genome and tracks to display by altering the 'url' parameter. 'trackfile' is the url of a file describing the custom tracks (beds/bigwigs) to display
mergePDF <- function(output="merged.pdf", sourcefiles=c("source1.pdf","source2.pdf","source3.pdf"))
{
    # create the command string and call the command using system()
    # merging command from http://hints.macworld.com/article.php?story=2003083122212228
    command=paste("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite",paste("-sOutputFile",output, sep="="), paste(sourcefiles, collapse=" "),sep=" ")
    try(system(command))
}


#  Reference: http://www.biostars.org/post/show/6132/batch-viewing-of-ucsc-browser-graphic/
screenshotUCSC <- function(url, trackfile, chr, start, end, filename) {
        oldpen <- options("scipen")
        options(scipen=100)
        temp <- readLines(paste(url, "&hgt.customText=", trackfile, "&position=",chr,":",start,"-",end, sep=""))
        #cat(temp,"\n")
        pdfurl <- paste("http://genome.ucsc.edu/trash",gsub(".*trash","",gsub(".pdf.*","",temp[grep(".pdf", temp, fixed=TRUE)][1])), ".pdf", sep="")
        cat(pdfurl,"\n");
        options(scipen=oldpen)
        download.file(pdfurl, filename, mode="wb", quiet=TRUE)
}

table2html <- function(df, caption="output.html", htmlfile="output.html", sparkline_col=1:ncol(df), sparkline_option='')
{
    # print header
    cat("<html>\n",  file=htmlfile)
    cat("\t<head>\n",  file=htmlfile, append=T)
    cat("\t\t<title>", caption, "</title>\n", file=htmlfile, append=T)
    cat("\t\t<script type=\"text/javascript\" src=\"http://code.jquery.com/jquery-1.8.0.js\"></script>\n",  file=htmlfile, append=T)
    cat("\t\t<script type=\"text/javascript\" src=\"http://omnipotent.net/jquery.sparkline/2.0/jquery.sparkline.js\"></script>\n",  file=htmlfile, append=T)
    cat("\t\t<script type=\"text/javascript\" src=\"http://zlab.umassmed.edu/~dongx/mylib.js\"></script>\n",  file=htmlfile, append=T)
    cat("\t\t<link rel=\"stylesheet\" type=\"text/css\" href=\"http://zlab.umassmed.edu/~dongx/stdtheme.css\" />\n",  file=htmlfile, append=T)
    cat("\t\t<script type=\"text/javascript\">$(function() {$('.sparkline').sparkline('html',{width:'40px',", sparkline_option, "});});</script>\n",  file=htmlfile, append=T)
    cat("\t</head>\n",  file=htmlfile, append=T)

    # print body
    cat("<body>\n",  file=htmlfile, append=T)
    cat("\t<table class='reference'>\n\t\t<caption><font size=5>", caption, "</font></caption>\n", file=htmlfile, append=T)
    # table header
    cat("\t\t<tr>\n\t\t\t<th>", file=htmlfile, append=T)
    cat("ID", sapply(colnames(df),function(x) as.character(x)), "sparkline", sep="</th>\n\t\t\t<th>", file=htmlfile, append=T)
    cat("</th>\n\t\t</tr>\n",  file=htmlfile, append=T)

    for(i in 1:nrow(df))
    {
        cat("\t\t<tr>\n\t\t\t<td>",  file=htmlfile, append=T)
        cat(rownames(df)[i], sapply(df[i,],function(x) as.character(x)), paste("<span class=\"sparkline\">",paste(df[i,sparkline_col],collapse=","),"</span>",collapse=""),
            sep="</td>\n\t\t\t<td>", file=htmlfile, append=T)
        cat("</td>\n\t\t</tr>\n",  file=htmlfile, append=T)
    }

    cat("\t\t</caption>\n\t</table>\n", file=htmlfile, append=T)

    # header/body end
    cat("</body>\n",  file=htmlfile, append=T)
    cat("</html>\n",  file=htmlfile, append=T)

}


getAGGREGATION2 <- function(plus, minus, genes,  TRIM=0, meanmedian='mean')
{
                          histdata=read.table(plus)
                          histdata=data.frame(histdata[,-1], row.names=histdata[,1])

                          p=histdata[rownames(genes)[genes$strand=='+'], ]
                          m=histdata[rownames(genes)[genes$strand=='-'], ]*-1

                          histdata=read.table(minus)
                          histdata=data.frame(histdata[,-1], row.names=histdata[,1])

                          m0=rbind(m,histdata[rownames(genes)[genes$strand=='+'], ])
                          p0=rbind(p,histdata[rownames(genes)[genes$strand=='-'], ]*-1)

                          # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
                          m0=m0[!is.na(rowSums(m0)),] # 79608 out of 79652 remained (for protein-cding gene)
                          p0=p0[!is.na(rowSums(p0)),]


                          # all genes
                          m1=m0; p1=p0;
                          # remove gene body bin, and 1% outlier
                          n1=nrow(m1)
                          if(meanmedian=="mean") {m1=apply(m1, 2, function(x) mean(x, trim=TRIM)); p1=apply(p1, 2, function(x) mean(x, trim=TRIM));}
                          if(meanmedian=="median") {m1=apply(m1, 2, median); p1=apply(p1, 2, median);}
                          #m1[41]=NA; p1[41]=NA
                          return(list(p1,m1))
}
# temp use for ENCODE.CSHL
getAGGREGATION2.cshl <- function(plus, minus, genes,  TRIM=0, meanmedian='mean')
{
                          histdata=read.table(plus)
                          histdata=data.frame(histdata[,-1], row.names=histdata[,1])

                          p=histdata[rownames(genes)[genes$strand=='+'], ]
                          m=histdata[rownames(genes)[genes$strand=='-'], ]*-1

                          histdata=read.table(minus)
                          histdata=data.frame(histdata[,-1], row.names=histdata[,1])

                          m0=rbind(m,histdata[rownames(genes)[genes$strand=='+'], ]*-1)
                          p0=rbind(p,histdata[rownames(genes)[genes$strand=='-'], ])

                          # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
                          m0=m0[!is.na(rowSums(m0)),] # 79608 out of 79652 remained (for protein-cding gene)
                          p0=p0[!is.na(rowSums(p0)),]


                          # all genes
                          m1=m0; p1=p0;
                          # remove gene body bin, and 1% outlier
                          n1=nrow(m1)
                          if(meanmedian=="mean") {m1=apply(m1, 2, function(x) mean(x, trim=TRIM)); p1=apply(p1, 2, function(x) mean(x, trim=TRIM));}
                          if(meanmedian=="median") {m1=apply(m1, 2, median); p1=apply(p1, 2, median);}
                          #m1[41]=NA; p1[41]=NA
                          return(list(p1,m1))
}

# aggregation plot for both-strand data
# filling the aear under lines with color
draw.plot2 <- function(p,m,ylim=range(p,m,na.rm=T), main="", xax=F, showlegend=F, ylab="", yax=F)
{
             #plot
             par(mar=c(2,2,1,1), lwd=2/3) # move to global setting
             plot(p, ylim=ylim, xaxt='n',yaxt='n', main='', type='n')
             if(xax==T) axis(1, at=c(1,20,61,81), labels=c('-1k','Start','End','+1k'), tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1.5)
             if(yax==T) axis(4, las=1, tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1.5)
             mtext(ylab, side=2, line=1, cex=2)
             mtext(main, side=3, line=0, cex=2, padj=0)

             #lines(c(1,81),c(0,0))
             abline(v=c(20,61), lty=2, col='darkgray')
             polygon(c(1,1:81,81),c(0,p[1:81],0), col='blue', border=NA)
             polygon(c(1,1:81,81),c(0,m[1:81],0), col='red', border=NA)
             if(showlegend==T) legend("topleft",c('sense strand','antisense strand'),col=c('blue','red'), bty='n', cex=1.5)
}


# Credit to Altuna's blog: http://zvfak.blogspot.com/2011/02/calling-bedtools-from-r.html
bedTools.2in<-function(functionstring="intersectBed",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

## Usage: pairs(USJudgeRatings[,c(2:3,6,1,7)], lower.panel=panel.smooth, upper.panel=panel.cor)
## Reference: http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=137
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))

    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex=cex, col=2)
}