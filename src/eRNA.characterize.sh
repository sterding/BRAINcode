## script to characterize eRNAs with the following features
#
#Features	            Data type	    Description
#========               =========       ==================================
#dis2TSS	            integer	        distance between the middle points of HiTNE and the nearest TSS, in bp. If HiTNE is located intronic, it's minus; otherwise it's positive. 
#RPKM	                float	        normalized expression level, calculated in the same way as RPKM.
#RPM	                float	        reads density at the summit position, normalized to total mapped reads in million. 
#readsCount	            integer	        raw reads count mapped to the HiTNE
#normCpG	            float	        normalized CpG score
#nTFBS	                integer	        number of distinct TFs bound to the HiTNE, based on ENCODE ChIPseq cluster (wgEncodeRegTfbsClusteredV3)
#P300	                boolean	        if the P300 binding site found in the HiTNE
#enhancer_CAGE	        boolean	        if overlap with any CAGE-defined permissive enhancers
#enhancer_histone	    boolean	        if overlap with any histone marks-defined enhancers (chromHMM states of E6|E7|E12 from substantial nigro)
#enhancer_VISTA	        boolean	        if overlap with any tested enhancers (positive enhancers from VISTA enhancer database)
#DNase	                float	        maximal DNase density of fetal brain from Roadmap Epigenomics
#bDNase	                boolean	        if overlap with DNase cluster from ENCODE (wgEncodeRegDnaseClustered V2)
#conservation	        float	        mean phastCons score for the HiTNE region
#bConserved2zf	        boolean	        if conserved to zebrafish (defined by existence of chain alignment between human and zebrafish)
#bHCNE	                boolean	        if overlapping with any HCNEs (HCNE_hg19_danRer7_70pc_50col from Ancora)
#GWAS	                integer	        number of GWAS SNPs in HiTNE
#bGWAS	                boolean	        if any GWAS SNPs in HiTNE
#eSNP	                integer	        number of eQTL SNPs in HiTNE
#beSNP	                boolean	        if any eQTL SNPs in HiTNE
#nHostgene	            integer	        number of HiTNEs in the host gene. 0 for intergenic HiTNE
#lenHostgene	        integer	        length of host gene. 0 for intergenic HiTNE
#enhancer_CAGE.HeLa	    boolean	        if overlap with any CAGE-defined enhancers in HeLa cell line  -- for in vitro validation purpose only
#enhancer_histone.HeLa	boolean	        if overlap with any histone marks-defined enhancers (chromHMM states of E6|E7|E12) from HeLa Cell -- for in vitro validation purpose only

input_bed="eRNA.bed"

ANNOTATION=$GENOME/Annotation/Genes

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

## dis2TSS
# ==================
fgrep -w gene $ANNOTATION/gencode.v19.annotation.gtf | sed 's/[;"]//g'  | awk '{OFS="\t"; print $1, $4-1, $5, $18"___"$10"___"$14, 0, $7}' > $ANNOTATION/gencode.v19.annotation.gtf.genes.bed
# have to deal with intronic and intergenic separately
# intronic ones (if located in two genes' intron, just randomly pick the first hit in the file.)
cat $ANNOTATION/gencode.v19.annotation.gtf.genes.bed | closestBed -a eRNA.bed -b - -d -t first | awk '{OFS="\t"; if($11==0) {tss=($10=="+")?$6:$7; d1=tss-$2; if(d1<0) d1=-d1; d2=tss-$3; if(d2<0) d2=-d2; print $4, (d1<d2)?-d1:-d2;}}' > /tmp/eRNA.tmp
# intergenic ones
cat $ANNOTATION/gencode.v19.annotation.gtf.genes.bed | closestBed -a eRNA.bed -b - -d -t first | awk '$11!=0' | cut -f1-4 | sort -k1,1 -k2,2n -u | closestBed -a - -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $ANNOTATION/gencode.v19.annotation.gtf.genes.bed) -D b -t first | awk '{OFS="\t"; if($11<0 && ($9+$11)<0) $11=-$11; print $4, $11}' >> /tmp/eRNA.tmp

sort -k1,1 /tmp/eRNA.tmp > eRNA.f1.dis2TSS.txt


### RPKM (Note: this RPKM is different from the normal RPKM, they might have a factor of read length difference)
# =================
# mean0: average over bases with non-covered bases counting as zeroes
#--------
inputBG=/data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bedGraph
bigWigAverageOverBed ${inputBG/bedGraph/bw} eRNA.bed stdout | cut -f1,5 | sort -k1,1 | awk '{OFS="\t"; print $1, $2*1000+0}' > eRNA.f2.RPKM.txt

### RPM 
# =================
bigWigAverageOverBed ${inputBG/bedGraph/bw} eRNA.bed stdout -minMax | cut -f1,8 | sort -k1,1  > eRNA.f3.RPM.txt

### raw reads count
# =================
# Note: use '-split' to exclude spliced reads when counting for intronic eRNA
for i in ~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.bam; do bsub -q normal -n 1 $HOME/neurogen/pipeline/RNAseq/modules/_get_readscount_per_region.sh $i eRNA.bed $i.eRNA.rawcount; done

## read in only the 80 subjects/samples (w/ genotype)
R
path="~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.bam.eRNA.rawcount"
IDs=read.table('~/neurogen/rnaseq_PD/results/merged/RNAseqID.wGenotyped.list',stringsAsFactors =F)[,1]
IDs=IDs[grep("^[HI].*_SNDA", IDs)]
EXP=data.frame(); id="locus"
for(i in Sys.glob(path)){
    ii=sub(".*run_output/(.*)/uniq.*","\\1", i);
    if(! (ii %in% IDs)) next;
    id=c(id, ii);    print(i)
    # read background
    expression=read.table(i, header=F) 
    if(ncol(EXP)==0) EXP=expression[,4:5] else EXP=cbind(EXP, expression[,5]); 
}
colnames(EXP)=id;
write.table(EXP, "eRNA.80samples.rawcount.xls", col.names=T, row.names=F, sep="\t", quote=F)
q('no')

# normalized CpG score
# ====================
bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed eRNA.bed -fo eRNA.seq.tab
getNormalizedCpGscore.awk eRNA.seq.tab | sort -k1,1 > eRNA.f5.CpG.tab

# overlap with VISTA enhancers
# ====================
feature=/PHShome/xd010/projects/PD/results/eRNA/externalData/VISTA/hg19.tested_regions.bed;
intersectBed -a eRNA.bed -b $feature -wao | awk '{OFS="\t"; print $4, ($6==-1)?"NONE":$8"___"$10"___"$11}' > eRNA.f11.VISTA.txt


# measure the splicing ratio & intron coverage of its host gene (if any), in order to remove the eRNAs likely being from the fuzzy genes (e.g. highly actively transcribed genes, co-transcriptional splicing? pre-mRNA? if intronic coverage > 30%)
# =================
# run premRNA.sh first
#--------
intersectBed -a eRNA.bed -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k15,15nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8-13,15 | sed 's/\t\./\t-1/g' > eRNA.premRNAratio

#intersectBed -a eRNA.bed -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wo | cut -f8 | sort | uniq -c | join -1 5 -2 2 -a 1 -e '0' -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.1' <(sort -k5,5 eRNA.premRNAratio) - | sed 's/ /\t/g' | sortBed > eRNA.premRNAratio.bed
#
#awk '{OFS="\t"; if($5!=-1) print $3-$2,$5,$11,$12;}' eRNA.premRNAratio.bed | sort -k2,2 | groupBy -g 2,3,4 -c 1 -o sum > genes.geneLength.eRNAcount.eRNAlength.tab

# how many eRNA are expected in the null distribution
sortBed -i /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/exons.meta.bed | mergeBed | cat - $GENOME/Annotation/Genes/hg19.gap.bed | bedtools shuffle -excl stdin -i eRNA.bed -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a - -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k15,15nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8-13,15 | sed 's/\t\./\t-1/g' > eRNA.shuffle.premRNAratio

## R
R
df1=read.table("eRNA.premRNAratio", header=F)
colnames(df1)=c("chr","start","end","ID", "host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length")
df2=read.table("eRNA.shuffle.premRNAratio", header=F)
colnames(df2)=c("chr","start","end","ID", "host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length")

library(ggplot2)
m <- ggplot(df1, aes(x = end-start))
m + geom_density() + scale_x_log10() +
geom_density(data=subset(df1,cov_intron>0.3), aes(x= end-start), color = "blue") +
geom_density(data=subset(df1,cov_intron<=0.3), aes(x= end-start), color = "red")
ggsave("eRNA.size.pdf", width=6, height=6)

require(plyr)
df=ddply(subset(cbind(df1, eRNAlength=df1$end-df1$start), host_gene!="-1"), c("host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length"), summarize, eRNAcount=length(hostGene_length), eRNAlength=sum(eRNAlength))
df=cbind(df, type="HiTNee")
df=rbind(df, cbind(ddply(subset(cbind(df2, eRNAlength=df2$end-df2$start), host_gene!="-1"), c("host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length"), summarize, eRNAcount=length(hostGene_length), eRNAlength=sum(eRNAlength)), type="shuffled"))

library(ggplot2)
ggplot(df, aes(eRNAcount, fill = type)) + geom_density(alpha = 0.2,adjust=3) + scale_x_log10() + theme_bw()+ theme(legend.position=c(1, 1), legend.justification=c(1,1))
ggsave("eRNAcount.in.hostgene.pdf", width=6, height=4)

df0=subset(df, host_gene %in% subset(df, type=='HiTNee' & eRNAcount>100)$host_gene)
df0$host_gene=gsub()
df0$host_gene2 = reorder(df0$host_gene, -df0$eRNAcount)
ggplot(df0, aes(x=host_gene2, y=eRNAcount, fill=type)) +
geom_bar(stat = "identity", position="dodge") +
coord_flip()+
theme_bw()+
theme(legend.position=c(1, 1), legend.justification=c(1,1), axis.text.y = element_text(size=5))

barplot(eRNAcount~type, df=subset(df, host_gene %in% subset(df, type=='HiTNee' & eRNAcount>100)$host_gene))
, aes(eRNAcount, fill=type)) + geom_bar() + theme_bw()+ theme(legend.position=c(1, 1), legend.justification=c(1,1))

df=df[with(df,order(type, -eRNAcount)),]
df=read.table("genes.eRNAcount.eRNAlength.tab")
colnames(df)=c("geneID","eRNAcount","geneLength")

pdf("eRNA.hostgene.pdf", width=12, height=8)
plot(eRNAcount~geneLength, df, col=ifelse(grepl("protein_coding",geneID),'#00000066','#ff000066'), pch=20)
text(df$geneLength, df$eRNAcount, labels=ifelse(df$eRNAcount>100,gsub("___ENSG.*","",df$geneID),""), cex=0.5, pos=4, offset=0.2, col=ifelse(grepl("protein_coding",df$geneID),'#000000','#ff0000'))

df=subset(df, eRNAcount>100)
df=df[with(df, order(-eRNAcount)),]
barplot(df$eRNAcount,
    names.arg =gsub("___ENSG.*","",df$geneID),
    ylab="eRNAs count in host gene",
    las=2, cex.names=.5,
    col=ifelse(grepl("protein_coding",df$geneID),"gray","red"),
    space=0.8,
    legend.text = c("protein-coding", "non-protein-coding"),
    args.legend = list(x = "topright", bty = 'n', fill=c("gray","red")))

dev.off()

q('no')


#############################################################
# eQTL of eRNA 
#############################################################
mkdir eQTL; cd eQTL
awk 'BEGIN{OFS="\t"; print "locus\tchr\ts1\ts2"}{print $4,$1,$2,$3}' ../eRNA.bed > eRNA.loci.txt
ln -fs ../eRNA.80samples.RPKM.xls eRNA.80samples.RPKM.xls
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_factor_analysis.R


#############################################################
# how much RNA-defined eRNA overlap with other enhancers --> venn diagram
#############################################################
# Roadmap Epigenomics enhancers (https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap)
# core chromatin states for substantial nigro: curl -s http://www.broadinstitute.org/~anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/E074_15_coreMarks_segments.bed > ../Segment/E074_15_coreMarks_segments.bed
cat ../Segment/E074_15_coreMarks_segments.bed | awk '$4~/E6|E7|E12/' | intersectBed -a eRNA.bed -b stdin -c | sort -k4,4  | cut -f4- > eRNA.overlap.txt

# CAGE-defined enhancers
intersectBed -a eRNA.bed -b ../CAGE/permissive_enhancers.bed -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

# DNase cluster  ( DNase cluster in V2: minimal score 500, at least 5 cell lines): http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeRegDnaseClustered
zcat ../DNase/wgEncodeRegDnaseClusteredV2.bed.gz | awk '$5>=500 && $4>=5' | intersectBed -a eRNA.bed -b stdin -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

# TFBS
intersectBed -a eRNA.bed -b ../TFBS/TFBS.distal.bed -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

# Conservation
curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_hg19_danRer7_70pc_50col.bed.gz | gunzip | intersectBed -a eRNA.bed -b stdin -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

## blood enhancer defined by H3k4me1/2/3  --> TO REMOVE
#intersectBed -a eRNA.bed -b ../Histone/blood.enhancers.PMID25103404.tab -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
#mv tmp.list eRNA.overlap.txt

# e2QTL
cut -f2 eQTL/final.cis.eQTL.FDR.05.xls | grep chr  | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | join -1 1 -2 2 -a 1 -e 0 -o '0,1.2,1.3,1.4,1.5,1.6,2.1' eRNA.overlap.txt - | sed 's/ /\t/g' > tmp.list
mv tmp.list eRNA.overlap.txt

# gene-eQTL
fgrep -w -f <(cut -f2 ~/neurogen/rnaseq_PD/results/eQTL/genes80samples/final.cis.eQTL.FDR.05.xls | grep -v SNP | sort -u) /data/neurogen/genotyping_PDBrainMap/eQTLMatrix/All.Matrix.SNP.ID | awk '{OFS="\t"; print $2,$3-1,$3,$1}' | intersectBed -a eRNA.bed -b stdin -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

## GWAS
# merge the recent PD-GWAS and NHGRI GWAS catalog SNPs
awk -v q="'" 'BEGIN{FS="\t"; OFS="\t";}{print $1,$2,$3,$4,"Parkinson"q"s disease", "25064009", $15}' ../GWAS/PD.GWAS.allsignificantSNP.bed > PDgwas.NHGRIgwascatalog.bed
awk 'BEGIN{FS="\t";OFS="\t";}{if($2!="PUBMEDID" && $13!="") print "chr"$12,$13-1,$13,$22,$8,$2,$28}' $GENOME/Annotation/GWASCatalog/gwascatalog2015Jan.txt >> PDgwas.NHGRIgwascatalog.bed 
# optional TODO: filter only the brain specific GWAS SNPs, or just PD GWAS
cut -f1-3 PDgwas.NHGRIgwascatalog.bed | sort -u | intersectBed -a eRNA.bed -b stdin -c | sort -k4,4 | cut -f5 | paste eRNA.overlap.txt - > tmp.list
mv tmp.list eRNA.overlap.txt

# fisher exact test to see the signifiance

R
df=read.table("eRNA.overlap.txt")
rownames(df)=df[,1]; df=df[,-1]
barplot(table(rowSums(df>0)), log='y', axes=F)
axis(2, 0:5, labels=10**c(0:5))

## robust eRNAs: at least one line of evidences from CAGE|Histone|TFBS|HCNE|eQTL|GWAS
echo -e "ID\tHistone\tCAGE\tDNase\tTFBS\tHCNE\te2QTL\teQTL\tGWAS\tn_evidences\tRPKM\thost_gene\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > robust_eRNAs.info.tab
awk '{s=0; for(i=2;i<=NF;i++) if($i>0) s++; OFS="\t"; print $0,s}' eRNA.overlap.txt | paste - <(cut -f2 eRNA.RPKM) | paste - <(cut -f5- eRNA.premRNAratio) | awk '$10>=3 && ($2>0 || $3>0 || $5>0 || $6>0 || $8>0 || $9>0)' >> robust_eRNAs.info.tab

echo -e "ID\tHistone\tCAGE\tDNase\tTFBS\tHCNE\te2QTL\teQTL\tGWAS\tn_evidences\tRPKM\thost_gene\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > permissive_eRNAs.info.tab
awk '{s=0; for(i=2;i<=NF;i++) if($i>0) s++; OFS="\t"; print $0,s}' eRNA.overlap.txt | paste - <(cut -f2 eRNA.RPKM) | paste - <(cut -f5- eRNA.premRNAratio) | awk '$10>0' >> permissive_eRNAs.info.tab

## get UCSC screenshot for eRNAs
grep chr robust_eRNAs.info.tab | cut -f1 | awk '{OFS="\t"; split($1,a,"_");print a[1],a[2],a[3],$1}' | Rscript ~/neurogen/pipeline/RNAseq/modules/_UCSC2PDF.R stdin
Rscript ~/neurogen/pipeline/RNAseq/modules/_table2html.R robust_eRNAs.info.tab robust.html robust_eRNAs


grep chr permissive_eRNAs.info.tab | cut -f1,10 | awk '{OFS="\t"; split($1,a,"_");print a[1],a[2],a[3],$1,$2}' | intersectBed -a stdin -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wo | awk '{OFS="\t"; print $0,$8-$7;}' |sort -k4,4 -k16,16nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f4,5,9,16 > eRNA.hostgene.tab

R
df=read.table("eRNA.hostgene.tab")
colnames(df)=c('ID','n_evidence','hostGene','length')
library("plyr")
df2=ddply(df, .(hostGene, length), summarize, n_eRNA=length(ID),sum_evidence=sum(n_evidence>=3))
df=ddply(df2, .(n_eRNA), summarize, n_gene=length(hostGene), mean_length=mean(length), mean_evidence=mean(sum_evidence))
pdf("eRNA.hostgene.pdf", width=10, height=4)
plot(df$n_gene, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
plot(df$mean_length, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
plot(df$mean_evidence, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
dev.off()



#############################################################
# measure eRNA, binarily and continously, with other features (e.g. CAGE, DNase, histone etc.)
#############################################################
# CAGE+
~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../CAGE/CAGE.fwd.bigwig eRNA.bed 1 max | sort -k1,1 > eRNA.onOtherFeatures.txt

# CAGE-
~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../CAGE/CAGE.rev.bigwig eRNA.bed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
mv tmp.list eRNA.onOtherFeatures.txt

# H3k4me1, H3K4me3 H3K27ac K3K27me3 H3K36me3 H3K9ac
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K9ac; do
    echo $i;
    ~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../Histone/Histone.SN.$i.bigwig eRNA.bed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
    mv tmp.list eRNA.onOtherFeatures.txt
done

# DNase, TFBS, Conservation
for i in DNase TFBS Conservation; do
    echo $i;
    ~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../$i/$i.bigwig eRNA.bed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
    mv tmp.list eRNA.onOtherFeatures.txt
done

# add header
echo -e "locus\tCAGE.fwd\tCAGE.rev\tH3K4me1\tH3K4me3\tH3K27ac\tH3K27me3\tH3K36me3\tH3K9ac\tDNase\tTFBS\tConservation" > tmp.list
cat eRNA.onOtherFeatures.txt >> tmp.list
mv tmp.list eRNA.onOtherFeatures.txt

## clustering
R
df=read.table("eRNA.onOtherFeatures.txt", header=T)
rownames(df)=df[,1]; df=df[,-1]

cage=df[,c('CAGE.fwd','CAGE.rev')]
library(pheatmap) # IF NOT, install.packages('pheatmap')
library("RColorBrewer")

pheatmap(log(1+df),scale='column',
    show_rownames=F,
    cluster_cols=F,
    clustering_distance_rows="correlation",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    filename="eRNA.onOtherFeatures.png")
)
