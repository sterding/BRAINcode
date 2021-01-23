###########################################
# R script to get the manhantten plot for eQTL SNPs
# Usage: Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R cis_eQTL_output.txt
# Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls ncRNA 
# Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls mRNA
# Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls all
# Author: Xianjun Dong
# Version: 0.0
# Date: 2016-Apr-25
###########################################
args<-commandArgs(TRUE)

eqtl_file_name=args[1] # eqtl_file_name="final.cis.eQTL.d1e6.p1e-2.xls"
genetype=args[2] # genetype="mRNA", ncRNA, or all

# dubug
# eqtl_file_name="final.cis.eQTL.d1e6.p1e-2.xls"; genetype='mRNA'

# resolution of x-axis
RES=5000000 

if(is.na(genetype)) genetype='all'

annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F, header=F); 
colnames(annotation) = c('chr','start','end','EnsID','score','strand','symbol','type')

message("# read eQTL data ... ")

df=read.table(eqtl_file_name, header=T, stringsAsFactors =F, sep="\t", quote="\"")
df$chr=sub("_.*","",df$gene); 
if(grepl("ENSG000", df$gene[1])){
  if(genetype=='mRNA') annotation=subset(annotation, type=="protein_coding")
  if(genetype=='ncRNA') annotation=subset(annotation, type!="protein_coding")
  df=subset(df, gene %in% annotation$EnsID)
  df$chr=annotation$chr[match(df$gene, annotation$EnsID)]
  #df$gene=do.call(paste, c(annotation[match(df$gene, annotation$EnsID),c('chr','start','end')],sep="_"))
}

message("# converting SNP coordinate ... ")

# only show chr1-22
df=subset(df, chr %in% paste0("chr",c(1:22)))

#head(df)
rcMpPvals=subset(df, p.value<=0.001, select=c(chr, SNP.pos, p.value))
## optional: for plotting efficiency, reduce the size of dataset by removing duplicated dots (or dots close enough to each other)
rcMpPvals$SNP.pos=round(rcMpPvals$SNP.pos/RES)     # precision = 10000 bp
rcMpPvals$p.value=round(-log10(rcMpPvals$p.value),2) # precision = 0.01
dim(rcMpPvals)
rcMpPvals=unique(rcMpPvals)
dim(rcMpPvals)

chrMax=0;bpMaxVec=c();
for (i in c(1:22)){ 
  ndx <- which(rcMpPvals[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  bpMaxVec=c(bpMaxVec, chrMax);
  rcMpPvals[ndx, 2] <- rcMpPvals[ndx, 2] + chrMax;
  chrMax = max(rcMpPvals[ndx, 2]);
  #print(c(i, length(ndx), chrMax))
}

bpMidVec <-c(); labels=c();
for (i in c(1:22)){
  ndx <- which(rcMpPvals[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  posSub <- rcMpPvals[ndx, 2]
  bpMidVec <- c(bpMidVec, ((max(posSub) - min(posSub))/2) + min(posSub))
  #labels = c(labels, paste0('chr',i))
  labels = c(labels, i)
}

message("# read the assocaited genes for top eQTL ... ")

topeSNPs=read.table(paste(dirname(eqtl_file_name),"final.cis.eQTL.d1e6.p1e-6.topSNPassociatedgene.xls",sep="/"), stringsAsFactors =F, header=F); 
head(topeSNPs)
colnames(topeSNPs)=c('chr','start','end','snp_gene','log10pvalue','associatedgene_eQTL_SNP','gene_type')
topeSNPs$gene = do.call(rbind,strsplit(topeSNPs$snp_gene,"|", fixed=T))[,2]
# only the eQTL on df dataset
if(genetype!='all') topeSNPs=subset(topeSNPs, gene %in% annotation$EnsID)
# only eQTL on chr1-22
topeSNPs=subset(topeSNPs, chr %in% paste0("chr",1:22))
topeSNPs=subset(topeSNPs, select=c(chr, end, log10pvalue, associatedgene_eQTL_SNP))

topeSNPs=merge(aggregate(log10pvalue~chr+end, data=topeSNPs, FUN=max),topeSNPs)
#topeSNPs=aggregate(associatedgene_eQTL_SNP~chr+end+log10pvalue, data=topeSNPs, FUN=paste0)

topeSNPs$end=round(topeSNPs$end/RES)     # precision = 10000 bp
topeSNPs$log10pvalue=round(topeSNPs$log10pvalue,2) # precision = 0.01

topeSNPs=merge(aggregate(log10pvalue~chr+end, data=topeSNPs, FUN=max),topeSNPs)
topeSNPs=unique(topeSNPs)

topeSNPs=aggregate(associatedgene_eQTL_SNP~chr+end+log10pvalue, data=topeSNPs, FUN=paste, collapse=',')
#topeSNPs[with(topeSNPs, order(chr, end, log10pvalue)), ]

message("# converting top eSNP coordinate ... ")

for (i in c(1:22)){ 
  ndx <- which(topeSNPs[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  topeSNPs[ndx, 2] <- topeSNPs[ndx, 2] + bpMaxVec[i];
}

message("# plotting ... ")

#rcMpPvals$chr[rcMpPvals$chr=="chrX"]="chr23"

library(ggplot2)
library(grid)
p=ggplot() +
  geom_point(data=rcMpPvals, aes(x=SNP.pos, y=p.value, colour=as.factor(as.numeric(sub("chr","",chr)) %% 2)),size=1, alpha=1) +
  scale_y_continuous(breaks=seq(3,12,3), limits=c(3,12)) + 
  scale_color_manual(values=c("0"='#888888',"1"='#222222'), name = "") +  
  theme_bw(base_size=15) + 
  theme(aspect.ratio = 0.18, legend.position='bottom', legend.key.size = unit(5, "points"), legend.text = element_text(size = 8*0.3514598)) + 
  guides(col = guide_legend(nrow = 6)) + 
  scale_x_continuous(expand = c(0.01, 0.01), labels=labels, breaks=bpMidVec) +
  geom_hline(yintercept=-log10(max(df$p.value[df$FDR<=0.05])), linetype=2, col='red', lwd=.4) +
  ggtitle(paste(eqtl_file_name, genetype,sep=".")) + xlab('') + ylab('-log10(P)')

message("# add SNP with top eQTL score")
p = p+geom_point(aes(x=end, y=log10pvalue), size=4, shape=18,data=topeSNPs)
message("# add labels (take the one with most pvalue then left-most coordinate, per disease per hostgene per chr)")
# ref: http://stackoverflow.com/questions/6289538/aggregate-a-dataframe-on-a-given-column-and-display-another-column
p=p+geom_text(aes(x=end, y=log10pvalue,size=2, label = as.character(associatedgene_eQTL_SNP)), angle = 90,size=8*0.3514598, vjust = 0.5, hjust=0, nudge_y=0.5, data=topeSNPs)

p

# see answer in: https://stackoverflow.com/questions/9992275/ggplot2-pdf-import-in-adobe-illustrator-missing-font-adobepistd/21756600#21756600
ggsave(height=3.5,width=7,dpi=100, filename=sub("xls",paste0(genetype,".only.pdf"),eqtl_file_name), useDingbats=FALSE)

message("# DONE")