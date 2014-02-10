# draw differentially expression plot based on data from cufflinks/htseq-count
# Usage:
# awk '{OFS="\t"; if($4=="=" && ) {printf("%s\t%s\t%s\t%s\t", $1,$2,$3,$4); for(i=5;i<=28;i++) {split($i, a,"|"); printf("%f\t",a[4]);} print "";}}' cuffcmp.tracking > cuffcmp.tracking.RPKM
# Rscript DE.R $output_dir PD Ct $HOME/projects/PD/data/gencode.v14.annotation.bed15
# 
args<-commandArgs(TRUE)
run_output=args[1]
pattern_PATIENT=args[2]
pattern_CONTROL=args[3]
annotation_path=args[4]

# read annotation
trans=read.table(annotation_path, header=F)
genes=unique(trans[,13:15])
rownames(genes)=genes[,1]; genes=genes[,-1]; colnames(genes)=c('type','gene_symbol')


###########################
# DESeq (readcount)
###########################
countsTable=genes
# case
for(i in dir(run_output, pattern = pattern_PATIENT)){
    print(i);
    df=read.table(paste(run_output, i, paste(i, 'hgseqcount.by.gene.tab', sep="."), sep='/'), header=F)
    if(nrow(df)==0) next;
    rownames(df)=df[,1]
    df0=intersect(rownames(countsTable), rownames(df))
    countsTable=cbind(countsTable[df0,], df[df0, 2, drop=F])
    colnames(countsTable)[ncol(countsTable)]=paste('readsCount',i,sep="_");
}
# control
for(i in dir(run_output, pattern = pattern_CONTROL)){
    print(i);
    df=read.table(paste(run_output, i, paste(i, 'hgseqcount.by.gene.tab', sep="."), sep='/'), header=T)
    rownames(df)=df[,1]
    df0=intersect(rownames(countsTable), rownames(df))
    countsTable=cbind(countsTable[df0,], df[df0, 2, drop=F])
    colnames(countsTable)[ncol(countsTable)]=paste('readsCount',i,sep="_");
}
countsTable <- countsTable[, c(-1,-2)]

# filter out genes with 0 values in all 10 samples
countsTable <- countsTable[apply(countsTable, 1, sum)>10,]  #  36932/50476 (73%) genes retained.
#countsTable <- countsTable[apply(countsTable<30, 1, sum)!=10,]  #  19732/50476 (39%) genes retained.
# add pseudocount 1 to avoid log2(0) issue later
countsTable <- countsTable+1

library(DESeq)
conds <- ifelse(grepl(pattern_PATIENT,colnames(countsTable)), pattern_PATIENT, pattern_CONTROL)
cds <- newCountDataSet( countsTable[grep("ENS", rownames(countsTable)),], conds )

# Normalisation
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
cds <- estimateDispersions(cds)
str(fitInfo(cds))
res <- nbinomTest( cds, pattern_PATIENT, pattern_CONTROL)

# cluster samples
vsd <- getVarianceStabilizedData( cds)
dists <- dist( t( vsd ) )
pdf("DESeq.sampleCluster.pdf", width=8, height=8)
heatmap( as.matrix( dists ), symm=TRUE, margin=c(14, 14))
dev.off()

# plot the DE
pdf("DESeq.MAplot.pdf", width=10, height=8)
plotDE <- function(res){
                     plot(res$baseMean, res$log2FoldChange,
                          xlab="baseMean= Average(mean.PATIENT + mean.CONTROL)",
                          ylab="log2(FoldChange)",
                          main="MA plot for the contrast Patients vs. Control",
                          log="x", pch=20, cex=0.5, cex.lab=2, cex.main=2, cex.axis=1.5,
                          col = ifelse( res$padj < .1, "red", "darkgray" )
                          )
}
par(mar=c(5, 5, 4, 2))
plotDE( res )
legend("bottomright", c("genes significant at 5% FDR","the rest"), cex=2, col=c('red','darkgray'), pch=20,  bty="n")
dev.off()

total=cbind(genes[res$id,], res, countsTable[res$id,])
write.table(total[order(total$padj),], "DESeq.xls", quote =F, sep = "\t", na="", row.names =F)

# DE genes
resSig=res[res$padj<0.1 & !is.na(res$padj),]
resSig=resSig[order(resSig$padj),]

# ================ down-regulated genes
updownregulatedgene = rbind(cbind(Regulation="down-regulated", head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ], 30 )),
                             cbind(Regulation="up-regulated",  head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ], 30 )))
# which genes are they?
sigsortgenes = cbind(genes[match(intersect(updownregulatedgene$id, rownames(genes)), rownames(genes)),],
                     updownregulatedgene[match(intersect(updownregulatedgene$id, rownames(genes)), updownregulatedgene$id),-2])
write.table(sigsortgenes, "DESeq.down.up.regulated.genes.txt",  quote =F, sep = "\t", row.names =F)

# ============================ top ones
sigsortgenes = cbind(genes[match(intersect(resSig$id, rownames(genes)), rownames(genes)),],
                     resSig[match(intersect(resSig$id, rownames(genes)), resSig$id),])
write.table(sigsortgenes[abs(sigsortgenes$log2FoldChange)>1,], "DESeq.topgenes.pajd0.1.FC2.txt", quote =F, sep = "\t", row.names =F)


