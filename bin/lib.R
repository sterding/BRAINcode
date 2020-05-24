# ===========================================================================
# code to make html for each DE result (see ~/rnaseq/src/_DE.R)
# ===========================================================================
makeNewImages <- function(df,...){
  imagename <- c()
  imagename3 <- c()
  for (i in 1:nrow(df)){
    ensId <- rownames(df)[i]
    symbol <- df$symbol[i]
    d <- data.frame(samples=colnames(assay(ntd)), expression_ntd=assay(ntd)[ensId,], expression_raw=assay(dds)[ensId,], cell=colData(dds)$Cell_type, dx=colData(dds)$Diagnosis, status=colData(dds)$Status)
    if(!file.exists(file.path('report/figures',paste('plot', ensId, symbol, 'txt', sep = ".")))) {
      write.table(d,file.path('report/figures',paste('plot', ensId, symbol, 'txt', sep = ".")),sep="\t", quote =F, row.names=F, col.names = T)
    }
    # 2 groups
    imagename[i] <- paste('plot', ensId, symbol, 'png', sep = ".")
    if(!file.exists(file.path('report/figures',imagename[i]))) {
      p=ggplot(d, aes(x=dx, y=expression_ntd, group=cell)) + 
        geom_jitter(size=1.5, position = position_jitter(width=.15)) +
        facet_wrap(~ cell, nrow=1) +
        stat_summary(fun.y=mean, geom="line", colour="red", size=0.8) +
        xlab("Diagnosis") + ylab("log2(counts+1)") + ggtitle(paste(ensId, "(",symbol,")"))
      png(filename = file.path('./report/figures', imagename[i]), height = 250, width = 600)
      print(p)
      dev.off()
    }
    # 3 groups
    imagename3[i] <- paste('plot3', ensId, symbol, 'png', sep = ".")
    if(!file.exists(file.path('report/figures',imagename3[i]))) {
      p=ggplot(d, aes(x=status, y=expression_ntd, group=cell)) + 
        geom_jitter(size=1.5, position = position_jitter(width=.15)) +
        facet_wrap(~ cell, nrow=1) +
        stat_summary(fun.y=mean, geom="line", colour="red", size=0.8) +
        xlab("Status") + ylab("log2(counts+1)") + ggtitle(paste(ensId, "(",symbol,")"))
      png(filename = file.path('./report/figures', imagename3[i]), height = 250, width = 600)
      print(p)
      dev.off()
    }
  }
  # df$Boxplot <- hwriteImage(paste0('figures/', imagename), 
  #                           link=paste0('figures/', imagename), 
  #                           table=FALSE, width=100)
  df$Boxplot <- hwrite("plot", link = paste0('figures/', imagename), table=F)
  df$Boxplot3 <- hwrite("plot3", link = paste0('figures/', imagename3), table=F)
  df$Rawdata <- hwrite("data", link = paste0('figures/', gsub("png","txt",imagename)), table=F)
  df$symbol <- hwrite(as.character(df$symbol), link = paste0("http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",as.character(rownames(df))), table=F)
  return(df)
}

# ===========================================================================
# GO enrichment analysis for a list of gene
# ===========================================================================
## Overlap problem in GO enrichment test: gene set are highly overlapped (e.g. Both gene sets, A and B, are enriched. However, B seems to be enriched only because of its overlap with A)
## topGO: offer solution to the overlap problem in gene set collection that are DAGs, by testing the leaf gene sets first and if a gene set is significant, remove its genes from its ancestor sets before testing the upper level sets.

# debug
# list_of_gene_names=DE_result$human_orth_symbol; allGenes=gene_names; topN=topN; pCutoff=Q_CUTOFF; type='all'; output=file.path(output_dir, input_DE_result)
topGOenrichment <- function(list_of_gene_names, allGenes, topN=100, pCutoff=0.001, type='all', output='input')
{
  if(!require("biomaRt")) { source("http://www.bioconductor.org/biocLite.R"); biocLite("biomaRt",ask=F);}
  if(!require("org.Hs.eg.db")) { source("http://www.bioconductor.org/biocLite.R"); biocLite("org.Hs.eg.db",ask=F);}
  if(!require("topGO")) { source("http://www.bioconductor.org/biocLite.R"); biocLite("topGO",ask=F);}
  if(!require("ggplot2")) install.packages('ggplot2',ask=F)
  
  # debug: 
  # list_of_gene_names=DE_result$human_orth_symbol;allGenes=gene_names
  
  input=toupper(unique(na.omit(list_of_gene_names))) # remove NA (for the otholog cases)
  all_genes =  factor(as.integer((allGenes %in% input)))  # input:1 rest: 0
  names(all_genes) <- allGenes
  table(all_genes)
  
  gt=data.frame();
  for(i in c("BP","CC", "MF")){
    GOdata <- new("topGOdata", ontology = i, 
                  allGenes = all_genes, 
                  description = "Test", 
                  annotationFun = annFUN.org, 
                  mapping = "org.Hs.eg.db", 
                  ID = "symbol")
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") # weight01 is a mixture between the elim and the weight algorithms
    # multiple test using all used GO (Ref: https://www.biostars.org/p/143083/) ---- not a practical sulution: nothing left after padj
    # gtable = GenTable(GOdata, classicFisher = resultFisher, topNodes = length(usedGO(GOdata)), numChar=100)
    gtable = GenTable(GOdata, classicFisher = resultFisher, topNodes = topN, numChar=100)
    gtable$classicFisher = as.numeric(gtable$classicFisher)
    gtable$classicFisher = ifelse(is.na(gtable$classicFisher), 1e-30, gtable$classicFisher);  # GenTable default trim anything <1e-30 as "<1e-30"
    gtable$padj = p.adjust(gtable$classicFisher, method='fdr')
    
    # add significant genes for each GO term
    allGO = genesInTerm(GOdata, whichGO = gtable$GO.ID) ##get all GOs and their genes from the topgo result.
    deGO = lapply(allGO,function(x) paste(x[x %in% input], collapse = ",")) ##get only DE genes in allGO
    gtable$Significant_genes = unlist(deGO) ## add to the gtable
    
    # sort by padj
    gtable = gtable[order(gtable$classicFisher),]
    ## filter with FDR and topN
    ## 8/8/18: only include gene set with size of [10,500]  (https://www.nature.com/articles/srep18871) (GSEA use [25,500], this paper use [10-200]: M.A. Mooney et al, 2015, Am J Med Genet B Neuropsychiatr Genet)
    gtable=head(subset(gtable, Annotated>=10 & Annotated<=500 & classicFisher <= pCutoff), topN)
    
    if(nrow(gtable)>0) gt=rbind(gt, data.frame(ontology=i, type=type, gtable))
  }    
  
  gt$Term <- factor(gt$Term, unique(as.character(gt$Term)))
  gt = gt[with(gt, order(ontology, -classicFisher)), ]
  
  write.table(gt, paste(output,type,"topGO.xls",sep="."), sep="\t", quote = F, col.names = NA, row.names = T)
  
  pdf(paste(output,type,"topGO.pdf",sep="."), width=9, height=nrow(gt)/5)
  p = ggplot(gt, aes(x = Term, y = -log10(classicFisher), fill=ontology)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_bw() + theme_classic() +
    ylab("GO terms") + xlab("-log10(Fisher's test P value)") + 
    ggtitle(paste("Top",topN,"enriched GO terms (p <",pCutoff,") in"), subtitle=paste(output,type))
  print(p)
  dev.off() 
  
}
# ===========================================================================
# ORA: over-representative analysis (GO, KEGG)
# ===========================================================================
ORA <- function(inputGenes, allGenes, topN=20, nCutoff=3, pCutoff=0.01, output='ORA.xls')
{
  inputGenes=toupper(unique(na.omit(inputGenes))) # remove NA (for the otholog cases)
  allGenes=toupper(unique(na.omit(allGenes))) # remove NA (for the otholog cases)
  
  gt=data.frame();
  #for(i in c("c5.bp","c5.cc", "c5.mf","c2.cp.kegg")){
  for(i in c("c2.cp.biocarta","c2.cp.reactome","c2.cp.kegg")){
    message(i);
    gmt.file = paste0("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/msigdb_v6.2/msigdb_v6.2_GMTs/",i,".v6.2.symbols.gmt")
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, utils::tail, -2)
    pathways = as.data.frame(cbind(do.call(rbind, lapply(pathwayLines,head,2)), geneset=pathways))
    fishtest <- function(x, y, all) {
      tb=table(X=factor(all %in% x[3][[1]], levels = c(TRUE, FALSE)), Y=factor(all %in% y, levels = c(TRUE, FALSE)))
      ft=fisher.test(tb)
      list(pvalue=ft$p.value, 
           OR=as.numeric(ft$estimate), 
           genelist=paste(intersect(x[3][[1]],y), collapse = ","), 
           n=length(intersect(x[3][[1]],y)))
    }
    oraRes=cbind(pathways[,1:2],do.call(rbind,apply(pathways, 1, fishtest, y=inputGenes, all=allGenes)))
    for(j in 1:ncol(oraRes)) oraRes[j] = unlist(oraRes[j])
    oraRes$FDR=p.adjust(oraRes$pvalue, method = 'fdr')
    oraRes$bonferroni=p.adjust(oraRes$pvalue, method = 'bonferroni')
    oraRes = oraRes[with(oraRes, order(pvalue)), ]
    if(nrow(oraRes)>0) gt=rbind(gt, data.frame(gene_set=i, oraRes))
  }
  
  # filter
  gt = subset(gt, FDR <= pCutoff & n > nCutoff)
  gt = gt[with(gt, order(gene_set, -OR)), ]

  write.table(gt, paste(output,"ORA.xls",sep="."), sep="\t", quote = F, col.names = NA, row.names = T)
  
  gt = gt %>% arrange(gene_set, -pvalue) %>% group_by(gene_set) %>% top_n(topN, -pvalue) %>% ungroup() 
  
  pdf(paste(output,"ORA.pdf",sep="."), width=9, height=nrow(gt)/5)
  p = mutate(gt,Term=factor(V1, unique(as.character(V1)))) %>%  
    ggplot(aes(x = Term, y = -log10(pvalue), fill=gene_set)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_bw() + theme_classic() +
    ylab("GO terms") + xlab("-log10(Fisher's test P value)") + 
    ggtitle(paste("Top",topN,"enriched GO terms (p <",pCutoff,") in"), subtitle=paste(output))
  print(p)
  dev.off() 
  
  return(gt)
}

# ===========================================================================
# SPIA analysis for a list of gene
# ===========================================================================
calcPathway <- function(file, bg, out="out", sp="hsa")
{
  
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

# ===========================================================================
# bigcor: Large correlation matrices in R
# http://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
# ===========================================================================
bigcor <- function(
  x, 
  fun = c("cor", "cov"), 
  size = 2000, 
  verbose = TRUE, 
  to.ffdf = FALSE, 
  ...)
{
  library(ff, quietly = TRUE)
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance"
  
  NCOL <- ncol(x)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
     
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))   
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)  
    
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  
  ## initiate time counter
  timeINIT <- proc.time() 
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat(STR, "of Block", COMB[1], "and Block", COMB[2], 
                     "(", length(G1), "x", length(G2), ") ... ")
    flush.console()
    RES <- FUN(x[, G1], x[, G2], ...)
    resMAT[G1, G2] <- RES
    resMAT[G2, G1] <- t(RES)   
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      cat(timeNOW[3], "s\n")
    }
    gc()
  }
  
  ## remove added colums from the beginning
  
   
  ## optionally convert to matrix format. Beware: takes long!
  if (to.ffdf) {
    if (verbose) cat("=> Converting to ffdf matrix ... ")
    resMAT <- as.ffdf(resMAT)
    timeNOW <- proc.time() - timeINIT
    if (verbose) cat(timeNOW[3], "s\n")
  }
  
  return(resMAT)
}

# a parallel version
bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
  library(ff, quietly = TRUE)
  require(doMC)
	if(ncore=="all"){
		ncore = multicore:::detectCores()
		registerDoMC(cores = ncore)
	} else{
		registerDoMC(cores = ncore)
	}
 
	NCOL <- ncol(x)
 
	## test if ncol(x) %% nblocks gives remainder 0
	if (NCOL %% nblocks != 0) {stop("Choose different 'nblocks' so that mod(ncol(x), nblocks) = 0!")}
 
	## preallocate square matrix of dimension
	## ncol(x) in 'ff' single format
	corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
 
	## split column numbers into 'nblocks' groups
	SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
 
	## create all unique combinations of blocks
	COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
	COMBS <- t(apply(COMBS, 1, sort))
	COMBS <- unique(COMBS)
 
	## iterate through each block combination, calculate correlation matrix
	## between blocks and store them in the preallocated matrix on both
	## symmetric sides of the diagonal
	results <- foreach(i = 1:nrow(COMBS)) %dopar% {
		COMB <- COMBS[i, ]
		G1 <- SPLIT[[COMB[1]]]
		G2 <- SPLIT[[COMB[2]]]
		if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
		flush.console()
		COR <- cor(MAT[, G1], MAT[, G2], ...)
		corMAT[G1, G2] <- COR
		corMAT[G2, G1] <- t(COR)
		COR <- NULL
	}
 
	gc()
	return(corMAT)
}



# cluster rows of X and graphically visualize it as a heatmap
plotClassification <- function(x)
{
    # x=df[,21:81]
    #d <- dist(x, method = "euclidean")
    d <- as.dist(1-cor(t(x)))
    
    hc <- hclust(d, method="average")
    plot(hc) # display dendogram
    image(d[[hc$order]])
    
    
    t<-read.table("histone/gencode_v3c_hg19_tr_with115_cshl_long_quantif.level12.protein_coding.TSS.wgEncodeBroadHistoneK562H3k79me2StdAln_2Reps.norm5.rawsignal.up2k.down10k.20bpBin",sep="\t")
    t.km<-kmeans(t,7)
    t.km.center<-apply(t.km$center,1,sum)
    write.table(t.km$cluster,"K562H3k79me2.cluster",quote=F,sep="\t",row.names=F,col.names=F)
    write.table(t.km.center,"K562H3k79me2.center",quote=F,sep="\t",row.names=F,col.names=F)
    gene.order<-order(order(order(t.km.center))[t.km$cluster])
    
    zmin<-0
    zmax<-8
    
    ## draw the 1174 and K562H3k79me2 as the figure2
    nheat<-1
    collist<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
    ColorRamp<-colorRampPalette(collist, bias=1)(10000)
    ColorLevels<-seq(to=zmax,from=zmin, length=10000)
    
    png("all.order.normSignal.K562H3k79me2.sub.png",width=6,height=8,units="in",res=300)
    par(mar=c(2,2,2,2))
    layout(matrix(seq(nheat+1),nrow=1,ncol=nheat+1),widths=c(rep(2,nheat),0.5),heights=rep(1,nheat+1))
    
    t<-as.matrix(t)
    t<-t*10000000/sum(t)
    t[t<zmin]<-zmin
    t[t>zmax]<-zmax
    ColorRamp_ex <- ColorRamp[round(1+(min(t)-zmin)*10000/(zmax-zmin)) : round( (max(t)-zmin)*10000/(zmax-zmin) )]
    image(1:ncol(t), 1:nrow(t),t(t[gene.order,]) , axes=FALSE, col=ColorRamp_ex,xlab="",ylab="")
    axis(2,at=c(1,cumsum(t.km$size[order(t.km.center)])),labels=F,tck=-0.02,mgp=c(3,0.2,0.2))
    
    image(1,ColorLevels,matrix(data=ColorLevels,ncol=length(ColorLevels),nrow=1),col=ColorRamp, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
    axis(4,seq(zmin,zmax,3),seq(zmin,zmax,3),mgp=c(mgp=c(3,0.2,0)))
    dev.off();
    
}

getAGGREGATION <- function(all, subtype,  TRIM=0, meanmedian='mean')
{
    subset=all;
    if(subtype!='all') subset=all[grep(subtype, rownames(all)),]

    subset=subset[!is.na(rowSums(subset)),] # remove trans with nan values from bigWigAverageBed (e.g. chr1 1 1)
    subset=subset[rowSums(subset)>0,] # Optional: remove genes with all zero (or unexpressed)

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
             plot(p, ylim=ylim, ylab='', xlab='', xaxt='n',yaxt='n', main='', type='S', lwd=2, col='red')
             abline(v=c(20,61), lty=2, col='darkgray')
             
             # old version
             #plot(p, ylim=ylim, ylab='', xlab='', xaxt='n',yaxt='n', main='', type='n')
             #lines(c(1,81),c(0,0))
             #abline(v=c(20,61), lty=2, col='darkgray')
             #polygon(c(1,1:81,81),c(0,p[1:81],0), col='red', border=NA)
             
             if(xax==T) axis(1, at=c(1,20,61,81), labels=c("-1k bp","5'","3'","+1k bp"), tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1.5)
             if(yax==T) axis(4, las=1, tck=.02, mgp=c(2,0.5,0), lwd=2/3, lwd.ticks=2/3,cex.axis=1)
             mtext(ylab, side=2, line=1, cex=1)
             mtext(main, side=3, line=1, cex=1, padj=0)
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