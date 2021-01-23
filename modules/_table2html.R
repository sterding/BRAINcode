# Script: converting tab to html
# Usage: Rscript _table2html.R eRNA.allinfo.tab index.html eRNA 1

args <- commandArgs(TRUE)
inputfile = args[1]
htmlfile = args[2]
caption = args[3]
toExcel = args[4]

zoomoutX=3

df=read.table(inputfile, header=T, stringsAsFactors =F)

# print header
cat("<html>\n",  file=htmlfile)
cat("\t<head>\n",  file=htmlfile, append=T)
cat("\t\t<title>", caption, "</title>\n", file=htmlfile, append=T)
cat("\t\t<script type=\"text/javascript\" src=\"http://code.jquery.com/jquery-1.8.0.js\"></script>\n",  file=htmlfile, append=T)
#cat("\t\t<script type=\"text/javascript\" src=\"http://omnipotent.net/jquery.sparkline/2.0/jquery.sparkline.js\"></script>\n",  file=htmlfile, append=T)
cat("\t\t<script type=\"text/javascript\" src=\"http://zlab.umassmed.edu/~dongx/mylib.js\"></script>\n",  file=htmlfile, append=T)
cat("\t\t<script type=\"text/javascript\" src=\"http://www.broadinstitute.org/~pouyak/encode-motif-disc/sorttable.js\"></script>\n",  file=htmlfile, append=T)
cat("\t\t<link rel=\"stylesheet\" type=\"text/css\" href=\"http://zlab.umassmed.edu/~dongx/stdtheme.css\" />\n",  file=htmlfile, append=T)
#cat("\t\t<script type=\"text/javascript\">$(function() {$('.sparkline').sparkline('html',{width:'40px',", sparkline_option, "});});</script>\n",  file=htmlfile, append=T)
cat("\t</head>\n",  file=htmlfile, append=T)

# print body
cat("<body>\n",  file=htmlfile, append=T)
cat("\t<table class='sortable reference'>\n\t\t<caption><font size=5>", caption, "</font></caption>\n", file=htmlfile, append=T)
# table header
cat("\t\t<tr>\n\t\t\t<th>", file=htmlfile, append=T)
cat(colnames(df), "ucsc screenshot", sep="</th>\n\t\t\t<th>", file=htmlfile, append=T)
cat("</th>\n\t\t</tr>\n",  file=htmlfile, append=T)

for(i in 1:nrow(df))
{
    message(i)
    cat("\t\t<tr>\n\t\t\t<td>",  file=htmlfile, append=T)
    pos=unlist(strsplit(as.character(df[i,1]),"_"))
    
    length=as.numeric(pos[3])-as.numeric(pos[2])
    chr=as.character(pos[1])
    start=ifelse(as.numeric(pos[2])>round(length*(zoomoutX-1)/2), as.numeric(pos[2])-round(length*(zoomoutX-1)/2), 1)
    end=as.numeric(pos[3])+round(length*(zoomoutX-1)/2)

    x=sapply(df[i,],function(x) as.character(x))
    genename = unlist(strsplit(x[12], "___"))[1]
    id=x[1]
    
    ## to Excel
    if(! is.null(toExcel)){
        df[i,1]=paste0("=hyperlink(\"http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sterding&hgS_otherUserSessionName=hg19_PD&db=hg19&position=",chr,":",start,"-",end,"\", \"", x[1],"\")")
        df[i,12]=paste0("=hyperlink(\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=",genename,"\", \"",x[12],"\")")
    }
        
    x[12]=ifelse(genename!="-1", paste0("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",genename," target='_blank'>",x[12],"</a>"), x[12])
    x[1] = paste0("<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sterding&hgS_otherUserSessionName=hg19_PD&db=hg19&position=",chr,":",start,"-",end," target='_blank'>",x[1],"</a>")
    
    cat(x, paste0("<a href=ucsc/eRNAucsc.",id,".png target='_blank'>png</a>"), sep="</td>\n\t\t\t<td>", file=htmlfile, append=T)
    cat("</td>\n\t\t</tr>\n",  file=htmlfile, append=T)
}

cat("\t\t</caption>\n\t</table>\n", file=htmlfile, append=T)

# header/body end
cat("</body>\n",  file=htmlfile, append=T)
cat("</html>\n",  file=htmlfile, append=T)

if(! is.null(toExcel)){
    write.table(df, paste0(htmlfile,".xls"), quote =F, sep ="\t", row.names =F)
}

