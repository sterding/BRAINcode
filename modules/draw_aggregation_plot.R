# script to draw aggregation plot of RNAseq signal for meta-gene
# input files: genes annotation file
# Usage: Rscript $HOME/neurogen/pipeline/RNAseq/modules/draw_aggregation_plot.R
# requirement: run bigWigAverageOverBed_81bins.sh before (see detail there)

source('/data/neurogen/pipeline/RNAseq/bin/lib.R')


trim=0.05
output_dir="/PHShome/xd010/neurogen/rnaseq_PD/results/aggregationPlot"
pattern="output_.*uniq.accepted_hits.normalized.bw"
gtype=c('all', 'protein_coding', 'lincRNA', 'pseudogene', 'processed_transcript', 'misc_RNA', 'miRNA', 'snoRNA', 'snRNA');

filenames = paste(output_dir, dir(output_dir, pattern=pattern), sep="/")  # for de novo piRNA genes

## merged categories (disease status + cell type, maybe also +batch+site)
mergedMP=list()

pdf("aggregationPlot_fixRange.pdf", height=2*length(gtype), width=6, colormodel='cmyk')
par(mfrow=c(length(gtype),1), oma=c(0,1,2,2))
for(fname in filenames){
    print(fname)
    genes=read.table(fname, header=F)
    genes=data.frame(genes[,-1], row.names=genes[,1])
    if(ncol(genes)!=81) next;
        
    tag=gsub(".*_of_(.*)_([:upper:]+).*_([:upper:]+)_([:digit:]).*", "\\1_\\2_\\3_\\4", fname)
    
    ## using individual range
    #for(gt in gtype){
    #    mp=getAGGREGATION(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')), trim);
    #    mergedMP[[tag]][[gt]]=colMeans(mergedMP[[tag]][[gt]], mp[[1]])
    #    draw.plot(mp[[1]], legend=paste(gt,"\n(n=",formatC(mp[[2]],format="d",big.mark=","),")",sep=""))
    #}

    # using the same range
    rangeY=c(0,0)
    MP=list()
    for(gt in gtype){
        mp=getAGGREGATION(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')), trim);
        mergedMP[[tag]][[gt]]=colMeans(mergedMP[[tag]][[gt]], mp[[1]])
        rangeY[1]=min(rangeY[1], range(mp[[1]], na.rm=T)[1])
        rangeY[2]=max(rangeY[2], range(mp[[1]], na.rm=T)[2])
        MP[[gt]]=mp
    }
    for(gt in gtype){
        draw.plot(MP[[gt]][[1]], ylim=rangeY, legend=paste(gt,"\n(n=",formatC(MP[[gt]][[2]],format="d",big.mark=","),")",sep=""))
    }
        
    mtext(gsub(".*_of_(.*)", "\\1", fname), outer = TRUE, cex = 1.2)
}
dev.off()

pdf("aggregationPlot_autoRange_merged.pdf", height=2*length(gtype), width=6, colormodel='cmyk')
par(mfrow=c(length(gtype),1), oma=c(0,1,2,2))
for(tag in names(mergedMP)){
    for(gt in gtype){
        mp=mergedMP[[tag]][[gt]]
        draw.plot(mp[[1]], legend=paste(gt,"\n(n=",formatC(mp[[2]],format="d",big.mark=","),")",sep=""))
    }
    mtext(tag, outer = TRUE, cex = 1.2)
}
dev.off()

quit("no")















## OLD VERSION [ DEPRECATED ]
QSUBOUPUT="../data"
filenames = paste(QSUBOUPUT, dir(QSUBOUPUT, pattern=pattern), sep="/")  # for de novo piRNA genes

png(paste("../result/aggregation.TSS.TTS",genetype,"png",sep="."), width=1200,height=1200)
par(mfrow=c(3,3))
for(type in c('CAGE','Degradome','PAS'))
{
             plus=filenames[grep(paste(type,'plus',sep=".*"), filenames)]
             minus=filenames[grep(paste(type,'minus',sep=".*"), filenames)]

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
             m1=apply(m1, 2, function(x) mean(x, trim=0.01)); m1[41]=NA
             p1=apply(p1, 2, function(x) mean(x, trim=0.01)); p1[41]=NA

             # HCP genes
             m2=m0[intersect(rownames(m0),hcp),]; p2=p0[intersect(rownames(m0),hcp),];
             # remove gene body bin, and 1% outlier
             n2=nrow(m2)
             m2=apply(m2, 2, function(x) mean(x, trim=0.01)); m2[41]=NA
             p2=apply(p2, 2, function(x) mean(x, trim=0.01)); p2[41]=NA

             # LCP genes
             m3=m0[intersect(rownames(m0),lcp),]; p3=p0[intersect(rownames(m0),lcp),];
             # remove gene body bin, and 1% outlier
             n3=nrow(m3)
             m3=apply(m3, 2, function(x) mean(x, trim=0.01)); m3[41]=NA
             p3=apply(p3, 2, function(x) mean(x, trim=0.01)); p3[41]=NA

             #plot
             draw.plot2(p1,m1,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (all genes; n=",formatC(n1,format="d",big.mark=","),")", sep=""))
             draw.plot2(p2,m2,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (HCP genes; n=",formatC(n2,format="d",big.mark=","),")", sep=""))
             draw.plot2(p3,m3,range(m1,m2,m3, p1,p2,p3, na.rm=T), paste(type," (LCP genes; n=",formatC(n3,format="d",big.mark=","),")", sep=""))

             print(type)
}
dev.off()