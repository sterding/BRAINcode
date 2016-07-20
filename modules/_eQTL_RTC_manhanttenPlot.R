###########################################
# R script to get the manhantten plot for eQTL SNPs
# Usage: Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R cis_eQTL_output.txt RTC.txt
# Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls ncRNA chr17:43583680-44506585
# Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls mRNA
# Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls all
# Author: Xianjun Dong
# Version: 0.0
# Date: 2016-Apr-25
###########################################
args<-commandArgs(TRUE)

eqtl_file_name=args[1] 
rtc_file_name=args[2] 
genetype=args[3] 

# debug
#eqtl_file_name="final.cis.eQTL.d1e6.p1e-2.xls"; rtc_file_name="final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls"; genetype='mRNA'; RES=5000000;


# resolution of x-axis
RES=5000000 

if(is.na(genetype)) genetype='all'

annotation=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", stringsAsFactors =F, header=F); 
colnames(annotation) = c('chr','start','end','EnsID','score','strand','symbol','type')

message("# read eQTL data ... ")

df=read.table(eqtl_file_name, header=T, stringsAsFactors =F, sep="\t", quote="\"")
df$chr=sub("_.*","",df$gene); 
if(grepl("ENSG000", df$gene[1])){
  if(genetype=='mRNA')  annotation=subset(annotation, type=="protein_coding")
  if(genetype=='ncRNA') annotation=subset(annotation, type!="protein_coding")
  df=subset(df, gene %in% annotation$EnsID)
  df$chr=annotation$chr[match(df$gene, annotation$EnsID)]
  #df$gene=do.call(paste, c(annotation[match(df$gene, annotation$EnsID),c('chr','start','end')],sep="_"))
}

message("# converting SNP coordinate ... ")

# only show chr1-22 (because LD, which RTC is based on, only has chr1-22)
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
  ndx <- which(rcMpPvals[, 1]==paste0('chr',i));
  if(length(ndx)==0) next;
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

message("# read RTC data ... ")

# highlight the ones with RTC>=0.85 and color by diseases
# setwd("~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA"); rtc_file_name="final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls"; 
rtc=read.table(rtc_file_name, header=F, sep="\t", quote ="",stringsAsFactors=F)
#head(rtc)
colnames(rtc)=c("eQTL_SNP", "assocaitedRNA", "pvalue_eQTL", "FDR_eQTL", "GWAS_SNP", "Complex_Trait", "RTC", "assocaitedRNA_hostgene", "assocaitedRNA_coord", "pos_GWAS_SNP", "pos_eQTL_SNP", "hostgene_GWAS_SNP", "hostgene_eQTL_SNP")
rtc$chr=sub("_.*","",rtc$assocaitedRNA_coord)
rtc$Complex_Trait=sub("\\(.*","",as.character(rtc$Complex_Trait))

## save simplied table
rtc0=unique(subset(rtc, RTC>=0.85, select=c(eQTL_SNP,chr,hostgene_eQTL_SNP,assocaitedRNA_hostgene,pvalue_eQTL,Complex_Trait,RTC)))
rtc0[is.na(rtc0)]="(intergenic)"
rtc0=merge(aggregate(pvalue_eQTL~eQTL_SNP+hostgene_eQTL_SNP+chr+assocaitedRNA_hostgene+Complex_Trait,data=rtc0,FUN=min),rtc0)
rtc0=subset(rtc0, select=c(eQTL_SNP,chr,hostgene_eQTL_SNP,assocaitedRNA_hostgene,pvalue_eQTL,Complex_Trait,RTC))
rtc0=rtc0[with(rtc0,order(chr,hostgene_eQTL_SNP,Complex_Trait, pvalue_eQTL)),]
write.table(rtc0,file=paste(rtc_file_name,"simplified.xls",sep="."),sep="\t", row.names = F, col.names = T, na="", quote=F)

# only those eQTL pairs in df (for cases that mRNA and ncRNA separately)
rtc=rtc[paste(rtc$chr,rtc$pos_eQTL_SNP,rtc$assocaitedRNA,sep="_") %in% paste(df$chr,df$SNP.pos,df$gene,sep="_"),]

rtc=subset(rtc, RTC>=0.85, select=c(chr,pos_eQTL_SNP,pvalue_eQTL,RTC,Complex_Trait,hostgene_eQTL_SNP))
rtc$Complex_Trait=factor(rtc$Complex_Trait, names(sort(table(rtc$Complex_Trait),decreasing =T)))
#levels(rtc$Complex_Trait) = paste0(names(table(as.character(rtc$Complex_Trait))), " (n=",as.numeric(table(as.character(rtc$Complex_Trait))),")")

rtc$pos_eQTL_SNP=round(rtc$pos_eQTL_SNP/RES)     # precision = 10000 bp
rtc$pvalue_eQTL=round(-log10(rtc$pvalue_eQTL),2) # precision = 0.01

library(made4) # source("http://bioconductor.org/biocLite.R"); biocLite("made4")
traits=getcol(length(levels(rtc$Complex_Trait)))
names(traits) = levels(rtc$Complex_Trait)

message("# converting RTC SNP coordinate ... ")

for (i in c(1:22)){ 
  ndx <- which(rtc[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  rtc[ndx, 2] <- rtc[ndx, 2] + bpMaxVec[i];
}

message("# plotting ... ")

#rcMpPvals$chr[rcMpPvals$chr=="chrX"]="chr23"

library(ggplot2)
library(grid)
p=ggplot() +
  geom_point(data=rcMpPvals, aes(x=SNP.pos, y=p.value, colour=as.factor(as.numeric(sub("chr","",chr)) %% 2)),size=1, alpha=1) +
  scale_y_continuous(breaks=seq(3,12,3), limits=c(3,12)) + 
  scale_color_manual(values=c("0"='#888888',"1"='#222222',traits), name = "") +  
  theme_bw(base_size=15) + 
  theme(legend.position='bottom', legend.key.size = unit(5, "points"), legend.text = element_text(size = 8)) + 
  guides(col = guide_legend(nrow = 6)) + 
  scale_x_continuous(expand = c(0.01, 0.01), labels=labels, breaks=bpMidVec) +
  geom_hline(yintercept=-log10(max(df$p.value[df$FDR<=0.05])), linetype=2, col='red', lwd=.4) +
  ggtitle(paste(eqtl_file_name, genetype,sep=".")) + xlab('') + ylab('-log10(P)')
if(nrow(rtc)>0){
  # add SNP with RTC score >= 0.85
  p = p+geom_point(aes(x=pos_eQTL_SNP, y=pvalue_eQTL, colour=as.factor(Complex_Trait)),size=4, shape=18,data=rtc)
  # add labels (take the one with most pvalue then left-most coordinate, per disease per hostgene per chr)
  # ref: http://stackoverflow.com/questions/6289538/aggregate-a-dataframe-on-a-given-column-and-display-another-column
  rtc=unique(subset(aggregate(pos_eQTL_SNP~chr+Complex_Trait+hostgene_eQTL_SNP+pvalue_eQTL,data=merge(aggregate(pvalue_eQTL~chr+Complex_Trait+hostgene_eQTL_SNP,data=rtc,FUN=max),rtc),FUN=min),select=c(pos_eQTL_SNP, pvalue_eQTL, hostgene_eQTL_SNP)))
  p=p+geom_text(aes(x=pos_eQTL_SNP, y=pvalue_eQTL,label = as.character(hostgene_eQTL_SNP)), size = 3, vjust = 0, nudge_y= 1, data=rtc)
}
p

# see answer in: https://stackoverflow.com/questions/9992275/ggplot2-pdf-import-in-adobe-illustrator-missing-font-adobepistd/21756600#21756600
ggsave(height=3.5,width=7,dpi=100, filename=sub("xls",paste0(genetype,".pdf"),eqtl_file_name), useDingbats=FALSE)

message("# DONE")