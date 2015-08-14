##############################################
## R script to test the feature of HiTNEs
## Author: Xianjun Dong
## Usage: Rscript $HOME/neurogen/pipeline/RNAseq/src/eRNA.characterize.R
## Date: Aug 10, 2015
## Version: 0.0
##############################################

features = read.table("~/eRNAseq/eRNA.characterize.xls", header=T, stringsAsFactors =F)
# add class
df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
features$class=ifelse(apply(df,1,sum)==0, 3, ifelse(df$f12.DNaseROADMAP==1, 1, 2))
features$gene_ENSID=sub(".*___(.*)___.*","\\1", features$f19.Hostgene)
# ===========================================================================
# GO enrichment for host genes of different HiTNE categories
# ===========================================================================

#source("http://www.bioconductor.org/biocLite.R"); biocLite(c("biomaRt", "topGO", "org.Hs.eg.db"))
library(biomaRt)
library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)

df <- read.table("/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F)
gene_names = df$V5[df$V6=="protein_coding"]

for(i in c("BP","CC", "MF")){
    
    all_genes =  factor(as.integer(gene_names %in% features$gene_ENSID))
    names(all_genes) <- sub("(.*)\\..*", "\\1", gene_names)
    GOdata <- new("topGOdata", ontology = i, allGenes = all_genes, geneSel = function(p) p < 0.01, description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Ensembl")
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    if(i=="BP") {
        gt=data.frame(ontology=i, type="all", GenTable(GOdata, classicFisher = resultFisher, topNodes = 10, numChar=100))
    } else {
        gt=rbind(gt, data.frame(ontology=i, type="all", GenTable(GOdata, classicFisher = resultFisher, topNodes = 10, numChar=100)))
    }

    for(k in unique(features$class))
    {
        all_genes =  factor(as.integer(gene_names %in% features$gene_ENSID[features$class==k]))
        names(all_genes) <- sub("(.*)\\..*", "\\1", gene_names)
        GOdata <- new("topGOdata", ontology = i, allGenes = all_genes, geneSel = function(p) p < 0.01, description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "Ensembl")
        resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
        gt=rbind(gt, data.frame(ontology=i, type=paste0("class",k), GenTable(GOdata, classicFisher = resultFisher, topNodes = 10, numChar=100)))
    }
}    

write.table(gt, "GOenrichment.genes.wHiTNE.xls", quote = F, col.names = NA, row.names = T)
gt=read.table("GOenrichment.genes.wHiTNE.xls", header=T)

gt$classicFisher = as.numeric(gt$classicFisher)
gt$classicFisher = ifelse(is.na(gt$classicFisher),1e-30,gt$classicFisher);  # GenTable default trim anything <1e-30 as "<1e-30"
gt=gt[with(gt, order(ontology,type,classicFisher)),]

pdf("GOenrichment.genes.wHiTNE.pdf", width=6, height=4, paper='us')
par(mar=c(4,16,2,2), mgp=c(2,.5,0));
for(i in c("BP","CC", "MF")){
    for(j in c('all','class1','class2','class3')){
        barplot(-log10(gt$classicFisher[gt$ontology==i & gt$type==j]), horiz=T, border=NA, names.arg=gt$Term[gt$ontology==i & gt$type==j], las=1, cex.names=.7, xlab="-log10(p-value)", main=paste0("(",i,",",j,")"))
    }
}
dev.off()  



#1. XYplot: Does the gene having more HITNEs tend to have more enhancers?

#2. Unique the HITNES and use pie charts to show the %
#3. Clustering the enhancers and see where the top clusters are located? Are they cooccuring with the neuronal genes?
#4. Genes with HITNEs vs. Genes without HITNEs, who has more enhancers? And who is longer?
#5. Mark the gene by HITNEs density, not just by HITNEs count. E.g. Any gene with 2x more HITNEs than expectation.
#6. Venn diagram of eHITNE (HITNEs that overlap with enhancers) and nHITNE (HITNEs that are located in the neuronal or synaptic genes)
#7. Any enrichment of  circRNAs in genes with enriched HITNEs?
#8. From the figure, It seems more enriched eRNA in intronic HITNEs than intergenic HITNEs.
#9. Add the eQTL and GWAS and disease gene info to select the "looks good" HITNEs


# length
df=read.table("eRNA.bed", header=F); 
df=data.frame(length=df$V3-df$V2, type="HiTNEs")
#df2=read.table("../CAGE/permissive_enhancers.bed", header=F, skip=1)
#df=rbind(df, data.frame(length=df2$V3-df2$V2, type="Enhancer RNAs (CAGE)"))
#
#df2=read.table("../CAGE/permissive_enhancers.bed", header=F, skip=1)
#df=rbind(df, data.frame(length=df2$V3-df2$V2, type="Enhancer RNAs (CAGE)"))
#df2=read.table("../CAGE/permissive_enhancers.bed", header=F, skip=1)
#df=rbind(df, data.frame(length=df2$V3-df2$V2, type="Enhancer RNAs (CAGE)"))
#df2=read.table("../CAGE/permissive_enhancers.bed", header=F, skip=1)
#df=rbind(df, data.frame(length=df2$V3-df2$V2, type="Enhancer RNAs (CAGE)"))

# dis2TSS : 94% of HiTNE are located in intragenic
pdf("eRNA.characterize.f01.dis2TSS.pdf", width=4, height=4, paper='usr')
d=hist(features$f01/1000, breaks=100, plot = F)
plot(d, col=ifelse(d$mids>0,'orangered','forestgreen'), border=NA, cex.axis=0.8, xlab="Distance to the nearest TSS (Kb)", ylab="Count", main="")
legend("topright", cex=0.8, paste0(c("intergenic", "intragenic")," (N=",format(table(features$f01>0), big.mark=','),"; ",round(100*table(features$f01>0)/nrow(features)),"%)"), col=c('forestgreen','orangered'), pch=15, bty='n')

# are they located in genes evenly? any 3' or 5' bias?
d=subset(features, f01.dis2TSS>0)
hist(100*d$f01/d$f21, breaks=80, xlab="Relative position to 5' TSS (%)", ylab="Count", border='orangered', main="", cex.axis=0.8, col='orangered')

d=density(100*d$f01/d$f21)
plot(d, xlab="Relative position to 5' TSS (%)", ylab="Density", main="", cex.axis=0.8, col='orangered')
polygon(d, col="orangered", border="orangered")

# intron length vs. HiTNE counts
#read random
rd=read.table("~/eRNAseq/intron.length.n_HITNE.random.txt", header=F, stringsAsFactors =F); colnames(rd)=c('intronID', 'len_intron','n_HITNE')
fd=read.table("~/eRNAseq/intron.length.n_HITNE.txt", header=F, stringsAsFactors =F);  colnames(fd)=c('intronID', 'len_intron','n_HITNE')
plot(fd$len_intron/1000, fd$n_HITNE, main="", xlim=range(rd$len_intron/1000, fd$len_intron/1000), col='orangered', cex.axis=0.8, cex=.6, pch=1, ylab="Number of HiTNEs in introns", xlab="Length of introns (Kb)")
points(rd$len_intron/1000, rd$n_HITNE, cex=.6, pch=3, col='darkblue')
legend("topleft", c("introns w/ HiTNEs","introns w/ randomly scrambled HiTNEs"), col=c('orangered','darkblue'), pch=c(1,3), cex=.8, bty='n')

#plot(fd$len_intron/1000, fd$n_HITNE, log='x', main="", xlim=range(rd$len_intron/1000, fd$len_intron/1000), col='orangered', cex.axis=0.8, cex=.6, pch=1, ylab="Number of HiTNEs in introns", xlab="Length of introns (Kb)")
#points(rd$len_intron/1000, rd$n_HITNE, cex=.6, pch=3, col='darkblue')
#legend("topleft", c("introns w/ HiTNEs","introns w/ randomly scrambled HiTNEs"), col=c('orangered','darkblue'), pch=c(1,3), cex=.8, bty='n')

rd=read.table("~/eRNAseq/Hostgene.length.nHITNErandom.metaExon.metaIntron.xls", header=F, stringsAsFactors =F); colnames(rd)=c('gene_symbol','gene_ENSID','gene_type', 'len_gene','n_HITNE','len_exon','len_intron')
fd=read.table("~/eRNAseq/Hostgene.length.nHITNE.metaExon.metaIntron.xls", header=F, stringsAsFactors =F);  colnames(fd)=c('gene_symbol','gene_ENSID','gene_type', 'len_gene','n_HITNE','len_exon','len_intron')

plot(fd$len_intron/1000, fd$n_HITNE, main="", xlim=range(rd$len_intron/1000, fd$len_intron/1000), col='orangered', cex.axis=0.8, cex=.6, pch=1, ylab="Number of HiTNEs in host gene", xlab="Total length of host gene introns (Kb)")
points(rd$len_intron/1000, rd$n_HITNE, cex=.6, pch=3, col='darkblue')
y = rd$n_HITNE; x <- rd$len_intron/1000
mod=lm(y ~ x)
newx <- seq(min(x), max(x), length.out=13377)
preds <- predict(mod, newdata = data.frame(x=newx), interval = 'confidence')
polygon(c(rev(newx), newx), as.numeric(c(rev(preds[ ,3]), preds[ ,2])), col = 'grey80', border = NA)
abline(mod, col='darkblue')
legend("topleft", c("genes w/ intronic HiTNEs","genes w/ randomly scrambled HiTNEs (y=0.0479*x+0.05)"), col=c('orangered','darkblue'), pch=c(1,3), cex=.8, bty='n')


## highlight the neuronal/synaptic genes
ng=read.table('neuronalgenes.Yunfei.txt', header=F, stringsAsFactors =F)
plot(fd$len_intron/1000, fd$n_HITNE, main="", xlim=range(rd$len_intron/1000, fd$len_intron/1000), col='orangered', cex.axis=0.8, cex=.6, pch=ifelse(fd$gene_symbol %in% ng$V1,19, 1), ylab="Number of HiTNEs in host gene", xlab="Total length of host gene introns (Kb)")
text(fd$len_intron/1000, fd$n_HITNE,labels=ifelse(fd$gene_symbol %in% ng$V1, fd$gene_symbol,""), cex=0.7, adj=c(0.5,1.5))

dev.off()

# ===========================================================================
# are they from PCR artifact?
# ===========================================================================
cpg=data.frame(region="HiTNEs", normalized.CpG=read.table("eRNA.f05.CpG.txt", stringsAsFactors =F)$V2, GC.content=read.table("eRNA.f05.CpG.txt", stringsAsFactors =F)$V3)
cpg=rbind(cpg, data.frame(region="promoter", normalized.CpG=read.table("promoters.f05.CpG.txt", stringsAsFactors =F)$V2, GC.content=read.table("promoters.f05.CpG.txt", stringsAsFactors =F)$V3))
cpg=rbind(cpg, data.frame(region="random", normalized.CpG=read.table("random.f05.CpG.txt", stringsAsFactors =F)$V2, GC.content=read.table("random.f05.CpG.txt", stringsAsFactors =F)$V3))

require('ggplot2')
ggplot(cpg, aes(normalized.CpG, colour = region)) + geom_density(adjust=5) +xlim(0,1)
ggplot(cpg, aes(GC.content, colour = region)) + geom_density(adjust=5) +xlim(0,1)



# ===========================================================================
# if the expression level of host gene is associated with the number of HITNEs
# ===========================================================================
expr=read.table("/PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.HCILB.uniq.xls", header=T, sep="\t", stringsAsFactors =F)
rownames(expr)=expr[,1]; expr=expr[,grep("FPKM",colnames(expr))]
expr=apply(expr, 1, mean, trim=.05)
# awk '{OFS="\t"; split($1,a,"_"); print a[1],a[2],a[3],$1,$2;}' eRNA.f02.RPKM.txt | intersectBed -a /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed -b - -wo | sort -k5,5 | groupBy -g 5 -c 12,12 -o count,mean > genes.nHiTNE.meanRPKM_HiTNE.txt
df=read.table("genes.nHiTNE.meanRPKM_HiTNE.txt"); rownames(df)=df[,1]; df=df[,-1]
fd=read.table("Hostgene.length.nHITNE.metaExon.metaIntron.xls", header=F, stringsAsFactors =F);  colnames(fd)=c('gene_symbol','gene_ENSID','gene_type', 'len_gene','n_HITNE','len_exon','len_intron')
fd=cbind(fd, medianFPKM=expr[fd$gene_ENSID], meanFPKM_HiTNE=as.numeric(as.character(df[fd$gene_ENSID,2])))

# grep ENSG gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed | cut -f4-  | sed 's/___/\t/g' | cut -f2,4- > /tmp/genes.rpkm_exon.rpkm_intron.cov_exon.cov_intron.splicing_ratio.txt
df=read.table("/tmp/genes.rpkm_exon.rpkm_intron.cov_exon.cov_intron.splicing_ratio.txt", sep="\t", stringsAsFactors =F); rownames(df)=df[,1];
fd=cbind(fd, rpkm_exon=df[fd$gene_ENSID,2], rpkm_intron=df[fd$gene_ENSID,3])

pdf("HITNE.vs.FPKM.pdf", paper='usr', width=4,height=4)
plot(log10(1e-3+fd$medianFPKM), fd$n_HITNE, xlab="mean FPKM of host gene", ylab='number of HiTNEs', pch=20, col='#FF4F0044')
plot(fd$medianFPKM, fd$meanFPKM_HiTNE, log='xy',pch=20, col='#FF4F0044', xlab='mean FPKM of host gene', ylab='mean FPKM of HiTNEs')
plot(fd$rpkm_exon+1e-3, 1e-3+fd$rpkm_intron, log='xy', asp=1, pch=20, col='#FF4F0044', xlab="RPKM of exons", ylab="RPKM of introns")
abline(a=0,b=1,lty=2)
dev.off()


# ===========================================================================
# P300 vs. TFBS
# ===========================================================================
pdf("eRNA.characterize.P300vsTFBS.pdf", width=10, height=4, paper='usr')
options(scipen=5)
barplot(table(features$f06.TFBS[features$f06>0]), log='y', ylim=c(.1,6000), cex.name=.6, border=NA, las=2, xlab="Number of TFBS per HiTNE", ylab="Count", yaxt="n")
barplot(table(features$f06.TFBS[features$f07>0]), log='y', ylim=c(.1,6000), cex.name=.6, border=NA, las=2, add=T, col='red', yaxt="n", xaxt='n')
axis(2, at=10**c(-1:3), labels=c(0,1,10,100,1000), las=2)
legend("topright", c("w/ P300 binding", "w/o P300 binding"), pch=15, col=c('red','gray'), bty='n')
dev.off()

