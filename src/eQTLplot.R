#awk '$1=="17" && $2<44506586 && $2>43583680' GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt | sort -k6,6 > snps.MAPTloci.pos
#cat *GTEx.csv | sed 's/"//g;s/,/\t/g;s/ //g' | cut -f2- | sort -k2,2 |join -1 2 -2 6 - <(sort -k6,6 snps.MAPTloci.pos) -o '1.1,0,1.3,1.4,2.2' | sed 's/ /\t/g' > eQTL.MAPTloci.txt

pdf('eQTLplot.MAPT.pdf', width=8.5, height=1.5); 


df=read.table("eQTL.MAPTloci.txt", header=F)
colnames(df) = c('gene','snp','pvalue','effectSize','pos')
df$gene = factor(df$gene, levels = c("MAPT","KANSL1"))
df$x=df$pos-43583680+1
require("ggplot2")
p <- ggplot(df, aes(x, -log10(pvalue)))
p=p + geom_point(aes(fill = gene), colour="#ffffff44", pch=21, size = 2) + xlim(1, 44506586-43583680+1) +  scale_x_continuous(expand = c(0,0))
#p + geom_point(aes(colour = factor(gene), alpha = effectSize))

print(p)

# eQTL from braincode
df=read.table("~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.FDRpt5.xls", header=F)
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
p = p + geom_point(aes(fill = gene), colour='#ffffff44', pch=21, size = 2) + scale_fill_discrete(drop = FALSE) + xlim(1, 44506586-43583680+1) +  scale_x_continuous(expand = c(0,0))

print(p)

dev.off()
