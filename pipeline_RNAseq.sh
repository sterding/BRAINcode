## A4: Pipeline for RNA-seq data analysis
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## date: 2010-12-05
## version: 1.0

######################################################
# Step0: setting
######################################################

inputdir="$HOME/projects/bu_neuro/data"
BOWTIE_INDEXES="/home/dongx/RDRIVE/hg19"

######################################################
# Step1: map reads to genome (Alignment)
######################################################
samlist=""
for src_fastqfile in `ls $inputdir/*.fastq`
do
    # Bowtie: allow 1 mismatch, soft-unique, 8-processors
    bowtie $BOWTIE_INDEXES/hg19 -v 1 -a --best --strata -m 1 -p 8 --quiet --sam ${src_fastqfile} > ${src_fastqfile%%.fastq}.sam

    # Tophat: To include reads that cross exon-exon boundaries
    #tophat -p 8 --max-multihits 2 --segment-mismatches 2 -o $src_fastqfile $BOWTIE_INDEXES/hg19 $src_fastqfile
    #mv $src_fastqfile/accepted_hits.bam ${src_fastqfile%%.fastq}.bam
    #samtools view ${src_fastqfile%%.fastq}.bam -o ${src_fastqfile%%.fastq}.sam
    #mv $src_fastqfile/junctions.bed ${src_fastqfile%%.fastq}.junctions.bed

    # Convert to sorted BAM format
    samtools import hg19.fa.fai ${src_fastqfile%%.fastq}.sam ${src_fastqfile%%.fastq}.bam
    #Sort everything
    samtools sort ${src_fastqfile%%.fastq}.bam ${src_fastqfile%%.fastq}_sorted
    #Delete temporary files
    mv ${src_fastqfile%%.fastq}_sorted.bam ${src_fastqfile%%.fastq}.bam
    #Create an index for fast searching
    samtools index ${src_fastqfile%%.fastq}.bam

    # Optional: generate BAM track files for UCSC browser
    [ -d $HOME/scratch/HD_tracks ] || mkdir -p $HOME/scratch/HD_tracks
    cp ${src_fastqfile%%.fastq}.bam $HOME/scratch/HD_tracks/
    cp ${src_fastqfile%%.fastq}.bam.bai $HOME/scratch/HD_tracks/
    # TODO: scp *.bam, *.bam.bai to zlab public track folder

    samlist="$samlist ${src_fastqfile%%.fastq}.sam"

    # Optional (only for Solution 2.2 below): using the htseq-count (http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)
    # to get reads count on each feature (e.g.gene, transcript). Below use gene
    htseq-count -m intersection-strict -t exon -i Parent -s no -q ${src_fastqfile%%.fastq}.sam $BOWTIE_INDEXES/hg19.gff > ${src_fastqfile%%.fastq}.hgseqcount.tab

    ## Optional: Run cufflinks for novel transcript discovery
    # run cufflinks to get FPKM
    cufflinks -o $src_fastqfile/cufflinks -j 0.05 -p 8 -Q 0 -M $BOWTIE_INDEXES/hg19.rRNA_MtRNA.gtf -r $BOWTIE_INDEXES/hg19.fa ${src_fastqfile%%.fastq}.sam
    # run cufflinks to get FPKM, only on hg19.gtf
    cufflinks -o $src_fastqfile/cufflinks_allknown -j 0.05 -p 8 -Q 0 -G $BOWTIE_INDEXES/hg19.gtf -M $BOWTIE_INDEXES/hg19.rRNA_MtRNA.gtf -r $BOWTIE_INDEXES/hg19.fa ${src_fastqfile%%.fastq}.sam

done

######################################################
# Step 2: map reads to transcriptome (Assemble)
######################################################
# Solution 2.1: ------ Using linux---------
cuffdiff -p 8 -o $HOME/scratch/HD_cuffdiff/ $BOWTIE_INDEXES/hg19.gtf $samlist

# R: deal with cuffdiff's output
fpkm = read.delim( "$HOME/scratch/HD_cuffdiff/genes.fpkm_tracking", header=T, stringsAsFactors=TRUE )
comp = read.delim( "$HOME/scratch/HD_cuffdiff/gene_exp.diff", header=T)
diffgenes = unique(comp[comp$significant=='yes', 1])

fpkm = fpkm[,grep("tracking|FPKM", colnames(fpkm))]
rownames( fpkm ) <- fpkm$tracking_id
fpkm <- fpkm[ , -1 ]

# only those genes with differentail expression between at least one case:control
fpkm=fpkm[diffgenes,]

# remove all-0 records
fpkm=fpkm[apply(fpkm, 1, sum)>0,]  # 12193/48065 = 25% left
fpkm[fpkm==0]=0.0001
Ct=fpkm[,grep("Ct", colnames(fpkm))]
HD=fpkm[,grep("HD", colnames(fpkm))]
ratio = log10(apply(HD,1,mean)/apply(Ct,1,mean))
hist(ratio[abs(ratio)>0.5], breaks=100)

# output genes with 100 times changes
topgenes = names(ratio)[abs(ratio) > 1]
fpkm[match(topgenes, fpkm$tracking_id),c(1,4)]
write.table(fpkm[match(topgenes, fpkm$tracking_id),c(1,4)], "../results/HD_cuffdiff/GENEID.TOP.txt")


# Solution 2.2: ------ Using htseqcount---------
# get reads count table
count_files=list.files(path="$inputdir", pattern="*.hgseqcount.tab$")
countsTable=c()
for(filename in count_files)
{
    tmp=read.delim(filename, header=F, stringsAsFactors=TRUE )
    countsTable=cbind(countsTable,tmp[,2])
    rownames( countsTable ) <- tmp[,1]
}
# colnames(countsTable) <- c("Ct_B6341", "HD_B4242", "HD_B3584", "Ct_B3732", "Ct_B6096", "HD_B4430", "HD_B4189")  # add label (optional)
conds <- substring(colnames(countsTable), 1, 2)

# R: using DESeq to get differentially expressed genes
library(DESeq)
cds <- newCountDataSet( countsTable[grep("ENST", rownames(countsTable)),], conds )

libsizes <- apply(countsTable, 2, sum)
sizeFactors(cds) <- libsizes
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "Ct", "HD")

sig=res[res$padj<0.1 & !is.na(res$padj),]

# top ones
sigsort=sig[with(sig, order(padj)), ]

# which genes are they?
geneanno = read.delim("/Users/dongx/projects/encode/data/Ensembl59.transcripts.annotation.tab")  # file download from Biomart for Ens v59
sigsortgenes = cbind(sigsort[match(intersect(sigsort$id, geneanno$Trans.ID), sigsort$id),],
                     geneanno[match(intersect(sigsort$id, geneanno$Trans.ID), geneanno$Trans.ID),c(1,2,8:11)])
write.table(head(sigsortgenes, n=200)$Gene.ID, "../results/DESeq/GENEID.TOP.txt")

# Solution 2.3 ------ Using R : Rsamtools---------

library(Rsamtools)
#Create a list of bam file object containing the short reads
bamlist=list()
src_files=list.files(path="$inputdir", pattern="*.bam$")
for(filename in src_files){
    #Since we do not know which strand the reads were originally transcribed,
    #so set the strand to be ambiguous
    tmp=readBamGappedAlignments(filename)
    bamlist[[length(bamlist)+1]] <- GRanges(seqnames=rname(tmp),
                                            ranges=IRanges(start=start(tmp),end=end(tmp)),
                                            strand=rep("*",length(tmp)))
}
names(bamlist)=src_files

# Ensembl annotation
library(GenomicFeatures)
txdb=makeTranscriptDbFromUCSC(genome="hg19",tablename="ensGene")

# summarize by gene and we choose to count all reads that fall within the body of the gene (including introns) as counting towards a genes count.
tx_by_gene=transcriptsBy(txdb,"gene")

#Initialize table of counts
toc=data.frame(rep(NA,length(tx_by_gene)))
for(i in 1:length(bamlist)){
    toc[,i]=countOverlaps(tx_by_gene,bamlist[[i]])
}
#Fix up labels
rownames(toc)=names(tx_by_gene)
colnames(toc)=names(bamlist)

######################################################
# Step 3: Differential Expression testing (Assess)
######################################################

# Normalization
library(edgeR)
norm_factors=calcNormFactors(as.matrix(toc))

# these scale factors are used as an offset in the negative binomial model.
DGE=DGEList(toc,lib.size=norm_factors*colSums(toc),group=gsub("[0-9].*","",colnames(toc)))

# Statistical Test
#  calculate a common dispersion parameter which represents the additional extra Poisson variability in the data.
disp=estimateCommonDisp(DGE)

tested=exactTest(disp)

######################################################
# Step 4: Gene Ontology testing (Annotation)
######################################################
library(goseq)
#Apply benjamini hochberg correction for multiple testing
#choose genes with a p-value less than .05
genes=as.integer(p.adjust(tested$table$p.value,method="BH") <.05)
names(genes)=row.names(tested$table)

#probability weighting function, which quantifies the length bias effect.
pwf=nullp(genes,"hg19","ensGene")

# Finally, we calculate the p-values for each GO category being over represented amongst DE genes.
GO.pvals=goseq(genes,pwf,"hg19","ensGene")

######################################################
# Step 5: Intersection with differentially enriched genes in Histone modification analysis (Annotation)
######################################################
# overlapped genes btw histone and RNA-seq
### ------ Hena's histone data
library(gdata)
his0 <- read.xls("../data/6HGT_vs_5nrm_agematch_062910.xls")
head(his0)  # look at the top few rows
his=his0[his0$DistanceToNearestTSS<100,]
his$log2_density_ratio=as.numeric(as.character(his$log2_density_ratio))  # de-factor

hisgenes = his$NearestTSS[abs(his$log2_density_ratio)>log2(2)]
write.table(hisgenes, "../results/DESeq/GENEID.his.txt", quote=F, col.names=F, row.names=F)

write.table(cbind(his[match(intersect(sigsortgenes$Gene.Name, hisgenes), his$NearestTSS),c(19,16)], sigsortgenes[match(intersect(sigsortgenes$Gene.Name, hisgenes), sigsortgenes$Gene.Name),]), sep="\t", "../results/DESeq/GENEID.overlapped.txt", quote=F, row.names=F)
