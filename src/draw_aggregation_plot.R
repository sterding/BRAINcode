# script to draw aggregation plot of RNAseq signal for meta-gene
# input files: genes annotation file
# Usage: Rscript $HOME/neurogen/pipeline/RNAseq/modules/draw_aggregation_plot.R
# requirement: run bigWigAverageOverBed_81bins.sh before (see detail there)

source('/data/neurogen/pipeline/RNAseq/bin/lib.R')

trim=0.05
output_dir="/PHShome/xd010/neurogen/rnaseq_PD/results/aggregationPlot"
pattern="output_.*uniq.accepted_hits.normalized.bw"
gtype=c('all', 'protein_coding', 'lincRNA', 'pseudogene', 'processed_transcript', 'housekeeping') #, 'misc_RNA', 'miRNA', 'snoRNA', 'snRNA');

# the longest Tx per gene
# from script as follow:
# sort -k4,4 /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.bed12 | sed 's/__ENST/\tENST/' | awk 'BEGIN{id=""; max=0;}{if($4!=id && id!="") {max=$6; print s; s=$0; } else if($6>max) {max=$6; s=$0;} id=$4; }END{print s;}' |  sed 's/\sENST/__ENST/'> /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.longestTx.bed12
longestTx=read.table("/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v19.annotation.longestTx.bed12")  # 57820 record

filenames = paste(output_dir, dir(output_dir, pattern=pattern), sep="/")  # for de novo piRNA genes

allMP=list()

## decay ratio
mRNAdecay=c();

mps=c();

pdf("aggregationPlot_fixRange.pdf", height=2*length(gtype), width=6, colormodel='cmyk')
par(mfrow=c(length(gtype),1), oma=c(0,1,2,2))
for(fname in filenames){
    print(fname)
    genes=read.table(fname, header=F)
    genes=data.frame(genes[,-1], row.names=genes[,1])
    if(ncol(genes)!=81) next;
        
    # only use the longest Tx as representative (in order to avoid the 3' end contamination)
    genes=genes[levels(longestTx[,4]),]
        
    tag=gsub(".*_of_([^.]*)..*", "\\1", fname)
    
    rangeY=c(0,0)
    MP=list()
    for(gt in gtype){
        ##make classification heatmap
        #plotClassification(genes, ifelse(gt=='all', gt, paste("__",gt,sep='')));
        
        # get aggregation data
        mp=getAGGREGATION(genes, ifelse(gt=='all', gt, ifelse(gt=="housekeeping", "^GAPDH_|^ACTB_|^KYNF_|^NEFL_|^B2M_", paste("__",gt,sep=''))), trim);
        allMP[[tag]][[gt]]=mp[[1]]
        mps=rbind(mps, mp[[1]]);
            
        # get decay ratio
        mRNAdecay=c(mRNAdecay, mean(mp[[1]][42:61])/mean(mp[[1]][21:40]))
                    
        rangeY[1]=min(rangeY[1], range(mp[[1]], na.rm=T)[1])
        rangeY[2]=max(rangeY[2], range(mp[[1]], na.rm=T)[2])
        MP[[gt]]=mp
    }
    for(gt in gtype){
        ## using the same range
        #draw.plot(MP[[gt]][[1]], ylim=rangeY, legend=paste(gt,"\n(n=",formatC(MP[[gt]][[2]],format="d",big.mark=","),")",sep=""))
        # using individual range
        draw.plot(MP[[gt]][[1]], legend=paste(gt,"\n(n=",formatC(MP[[gt]][[2]],format="d",big.mark=","),")",sep=""))
    }
        
    mtext(gsub(".*_of_(.*)", "\\1", fname), outer = TRUE, cex = 1.2)
}
dev.off()

mps=cbind(sampleID=rep(gsub(".*_of_([^.]*)..*", "\\1", filenames), each=length(gtype)), gType=rep(gtype, length(filenames)), mps)

mRNAdecay=matrix(mRNAdecay, ncol=length(gtype), byrow=T, dimnames=list(gsub(".*_of_([^.]*)..*", "\\1", filenames), gtype))

## decay ratio vs. RIN
# load covaraince table
covarianceTableURL="~/neurogen/rnaseq_PD/results/DE_DESeq2/covariances.tab"
covarianceTable=read.table(covarianceTableURL, header=T)
rownames(covarianceTable)=covarianceTable$sampleName;
covarianceTable=covarianceTable[,-1]
df=cbind(mRNAdecay, covarianceTable[match(rownames(mRNAdecay), rownames(covarianceTable)),])

df=df[grep("ILB_BN10-90_SNDA_4|ILB_BN99-54_SNDA_4", rownames(df), invert=T),] # remove two outliers

pdf("decayratio_vs_RIN.pdf", height=6, width=6, colormodel='cmyk')
plot(housekeeping~RIN, df, ylab="decay ratio", main="housekeeping genes"); legend("topleft", paste("Spearman's rho",round(cor(df$housekeeping, df$RIN, method='spearman'),3), sep="=")); abline(lm(housekeeping~RIN, df), col='red')
plot(protein_coding~RIN, df, ylab="decay ratio", main="protein_coding genes"); legend("topleft", paste("Spearman's rho",round(cor(df$protein_coding, df$RIN, method='spearman'),3), sep="=")); abline(lm(protein_coding~RIN, df), col='red')
plot(all~RIN, df, ylab="decay ratio", main="all genes"); legend("topleft", paste("Spearman's rho",round(cor(df$all, df$RIN, method='spearman'),3), sep="=")); abline(lm(all~RIN, df), col='red')
dev.off()

# save workspace into 
save.image("agg.RData")

# make the plot for merged group
## merged categories (disease status + cell type, maybe also +batch+site)
mergedMP=list()
for(sname in sort(names(allMP))){
    tag=gsub("(.*)_([[:upper:]]+).*_([[:upper:]]+)_([[:digit:]]).*", "\\1_\\2_\\3_\\4", sname)
    for(gt in gtype){
        if(is.null(mergedMP[[tag]][[gt]])) mergedMP[[tag]][[gt]]=allMP[[sname]][[gt]]
        else mergedMP[[tag]][[gt]]=colMeans(rbind(mergedMP[[tag]][[gt]], allMP[[sname]][[gt]]))
    }
}
pdf("aggregationPlot_autoRange_merged.pdf", height=2*length(gtype), width=6, colormodel='cmyk')
par(mfrow=c(length(gtype),1), oma=c(0,1,2,2))
for(tag in sort(names(mergedMP))){
    for(gt in gtype){
        draw.plot(mergedMP[[tag]][[gt]], legend=gt)
    }
    mtext(tag, outer = TRUE, cex = 1.2)
}
dev.off()


GENES=list();
for(fname in filenames){
    print(fname)
    genes=read.table(fname, header=F)
    genes=data.frame(genes[,-1], row.names=genes[,1])
    if(ncol(genes)!=81) next;
    GENES[[gsub(".*_of_(.*).uniq.*", "\\1", fname)]]=genes;
}

# make the plot for specific gene
queryGene="ENST00000394989" # e.g. SNCA genes

df=data.frame()        
for(n in names(GENES)){
    print(n)
    genes=GENES[[n]];
    d=genes[grep(queryGene,rownames(genes)),]
    rownames(d)=n #paste(n, gsub("(.*)__.*","\\1",rownames(d)), sep="___")
    df=rbind(df, d)
}

rownames(df)=gsub("(.*)___.*","\\1",rownames(df))

x=df[,21:61]
    #d <- dist(x, method = "euclidean")
    d <- as.dist(1-cor(t(x)))
    hc <- hclust(d, method="average")
    pdf(paste(queryGene, "tree.pdf", sep="."), width=5, height=15)
    par(mar=c(4,2,2,9)); plot(as.dendrogram(hc),horiz=T, cex=.4)
    #plot(hc, hang=-1)
    
    image(d[[hc$order]])
    
plotClassification(df[,21:])
# 