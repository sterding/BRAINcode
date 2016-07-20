# Rscript to make manhatten like plot for eQTL in MAPT locus, using GTEx, BRAINCODE result


### Gene eQTL from GTEx
## ===================

#cd ~/eRNAseq/externalData/GTEx/
#axel -a  -n 5 http://www.gtexportal.org/static/datasets/gtex_analysis_v6/single_tissue_eqtl_data/GTEx_Analysis_V6_eQTLs.tar.gz
#tar -zxvf GTEx_Analysis_V6_eQTLs.tar.gz
#awk '{OFS="\t"; if($14=="17" && $15<44506586 && $15>43583680) print $27,$23,$6,$3,$15}' Brain_Cerebellum_Analysis.snpgenes > eQTL.MAPTloci.all.txt
#awk '$1=="17" && $2<44506586 && $2>43583680' ~/eRNAseq/externalData/GTEx/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt | sort -k6,6 > snps.MAPTloci.pos
#cat ~/eRNAseq/externalData/GTEx/*GTEx.csv | sed 's/"//g;s/,/\t/g;s/ //g' | cut -f2- | sort -k2,2 |join -1 2 -2 6 - <(sort -k6,6 snps.MAPTloci.pos) -o '1.1,0,1.3,1.4,2.2' | sed 's/ /\t/g' > eQTL.MAPTloci.txt

df=read.table("~/eRNAseq/externalData/GTEx/eQTL.MAPTloci.all.txt", header=F)
colnames(df) = c('gene','snp','pvalue','effectSize','pos')
table(as.character(df$gene))

# ARL17B     CRHR1-IT1        KANSL1    KANSL1-AS1       LRRC37A     LRRC37A4P          MAPT      MAPT-AS1         NSFP1 
# 4             1             4          3555          3420          3650          3478          3305            56 
# RP11-259G18.1 RP11-259G18.2 RP11-259G18.3 RP11-707O23.5  RP11-798G7.5  RP11-798G7.8        SPPL2C 
# 3548          3590          3613          3568          3517          2795          2928 

df$gene = factor(df$gene, levels = c("MAPT","KANSL1"))
df$x=df$pos-43583680+1
require("ggplot2")
p <- ggplot(df, aes(x, -log10(pvalue)))
p + geom_point(aes(fill = gene), colour="#ffffff44", pch=21, size = 2) + xlim(1, 44506586-43583680+1) +  scale_x_continuous(expand = c(0,0))
#p + geom_point(aes(colour = factor(gene), alpha = effectSize))

### Gene eQTL from BrainCODE
## ===================
df=read.table("~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.bed6", header=F, stringsAsFactors = F)
df=read.table("~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.bed6", header=F, stringsAsFactors = F)
df=cbind(df, do.call(rbind,strsplit(as.character(df$V4),"\\|")))
colnames(df) = c('chr','SNP.start','SNP.end','eQTL_pair','p.value',	'FDR', 'snp','gene')
df$gene=as.character(df$gene); df$snp=as.character(df$snp)
genes=read.table("~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed", header=F, stringsAsFactors = F)
colnames(genes) = c('chr','start','end','geneID','score','strand','geneSymbol','type')
df=cbind(df, genes[match(df$gene, genes$geneID),c('start','end','geneSymbol')])
df$x=df$SNP.end-43583680+1
df=subset(df, chr=='chr17' & SNP.end<=44506586 & SNP.end>=43583680 & start<=44506586 & end>=43583680)
table(df$geneSymbol)
# FDR < 0.05
# LRRC37A4P     NSFP1 
# 2366        31 
# p< 0.01
# ARL17B        DND1P1        KANSL1    KANSL1-AS1       LRRC37A     LRRC37A4P          MAPT         NSFP1 RP11-105N13.4 
# 28          2461            38          2498            47          2520             1            40            48 
# RP11-293E1.1 RP11-669E14.4 RP11-707O23.5  RP11-798G7.5           STH 
# 1            20          2510            54             1 

# range of MAPT and KANSL1 chr17:43,969,657-44,313,597
df=subset(df, chr=='chr17' & start<=44313597 & end>=43969657 & SNP.end<=44506586 & SNP.end>=43583680)
df$gene=ifelse(df$start<=44105699 & df$end>=43971748, 'MAPT','KANSL1')
df=subset(df, select=c(gene, x, p.value))
# select the most significant p-value per gene per SNP, for a SNP may be linked to two different HTNEs from the same gen locus
library(dplyr)
df %>% group_by(gene, x) %>% summarise_each(funs(min(., na.rm=TRUE)), p.value)

df$gene = factor(df$gene, levels = c("MAPT","KANSL1"))

p <- ggplot(df, aes(x=x, y=jitter(-log10(p.value),0.2,0))) 
p + geom_point(aes(fill = gene), colour='#ffffff44', pch=21, size = 2) + scale_fill_discrete(drop = FALSE) + xlim(1, 44506586-43583680+1) +  scale_x_continuous(expand = c(0,0))


### HTNE eQTL from BrainCODE
## ===================
df=read.table("~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls", header=F)
df=cbind(df, do.call(rbind,strsplit(as.character(df$V2),"_")))
colnames(df) = c('snp','HTNE','beta',	't.stat',	'p.value',	'FDR',	'SNP.pos','chr','start','end')
df$start = as.numeric(as.character(df$start)); df$end = as.numeric(as.character(df$end))
df$x=df$SNP.pos-43583680+1
# range of MAPT and KANSL1 chr17:43,969,657-44,313,597
df=subset(df, chr=='chr17' & start<=44313597 & end>=43969657 & SNP.pos<=44506586 & SNP.pos>=43583680)
df$gene=ifelse(df$start<=44105699 & df$end>=43971748, 'MAPT','KANSL1')
df=subset(df, select=c(gene, x, p.value))
# select the most significant p-value per gene per SNP, for a SNP may be linked to two different HTNEs from the same gen locus
library(dplyr)
df %>% group_by(gene, x) %>% summarise_each(funs(min(., na.rm=TRUE)), p.value)

df$gene = factor(df$gene, levels = c("MAPT","KANSL1"))

p <- ggplot(df, aes(x=x, y=jitter(-log10(p.value),0.2,0))) 
p + geom_point(aes(fill = gene), colour='#ffffff44', pch=21, size = 2) + scale_fill_discrete(drop = FALSE) + xlim(1, 44506586-43583680+1) +  scale_x_continuous(expand = c(0,0))

ggsave(height=1.5,width=8.5,dpi=200, filename='eQTLplot.MAPT.pdf', useDingbats=FALSE)