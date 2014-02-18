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
        
    tag=gsub(".*_of_(.*)_([[:upper:]]+).*_([[:upper:]]+)_([[:digit:]]).*", "\\1_\\2_\\3_\\4", fname)
    
    ## using individual range
    #for(gt in gtype){
    #    mp=getAGGREGATION(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')), trim);
    #    if(is.null(mergedMP[[tag]][[gt]])) mergedMP[[tag]][[gt]]=mp[[1]]
    #    else mergedMP[[tag]][[gt]]=colMeans(rbind(mergedMP[[tag]][[gt]], mp[[1]]))
    #    draw.plot(mp[[1]], legend=paste(gt,"\n(n=",formatC(mp[[2]],format="d",big.mark=","),")",sep=""))
    #}

    # using the same range
    rangeY=c(0,0)
    MP=list()
    for(gt in gtype){
        #make classification heatmap
        plotClassification(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')));
        
        # get aggregation data
        mp=getAGGREGATION(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')), trim);
        if(is.null(mergedMP[[tag]][[gt]])) mergedMP[[tag]][[gt]]=mp[[1]]
        else mergedMP[[tag]][[gt]]=colMeans(rbind(mergedMP[[tag]][[gt]], mp[[1]]))
        
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

# save workspace into 
save.image("agg.RData")

# make the plot for merged group
pdf("aggregationPlot_autoRange_merged.pdf", height=2*length(gtype), width=6, colormodel='cmyk')
par(mfrow=c(length(gtype),1), oma=c(0,1,2,2))
for(tag in sort(names(mergedMP))){
    for(gt in gtype){
        draw.plot(mergedMP[[tag]][[gt]], legend=gt)
    }
    mtext(tag, outer = TRUE, cex = 1.2)
}
dev.off()

# 