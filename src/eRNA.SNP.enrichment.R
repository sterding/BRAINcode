# The Fisher exact test was performed to test the odds ratio of trait association with enhancers vs non-enhancers taking into account the genomic distribution of SNPs.  
# OR = (A/B) / (C/D), where
# 
# A: number of trait-associated SNPs in enhancers
# B: number of non-trait-associated SNPs in enhancers (from dbSNP - GWAS SNPs for that trait)
# C: number of trait-associated SNPs distal to enhancers
# D: number of non-trait-associated SNPs distal to enhancers (from dbSNP - GWAS SNPs for that trait)

args<-commandArgs(TRUE)
type=args[1]

require(ggplot2)

# setwd("~/eRNAseq/HCILB_SNDA"); type="PLINK"

s=read.table(paste0("SNP.",type,".counts.summary"), header=F,row.names=1); 
results=data.frame();
for(i in c('HTNE','pHTNE','exon','promoter','random')){
  n1=s[i,1]; n2=s['all',1];
  all=read.table(paste0("SNP.",type,".count.all")); rownames(all)=all[,1]
  x=read.table(paste0("SNP.",type,".count.",i)); rownames(x)=x[,1]
  df=cbind(x, all[rownames(x),2]); df=df[,-1]; colnames(df)=c('observed','all')
  results=rbind(results, cbind(Disease_or_Trait=rownames(df), 
                               df, 
                               pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), 
                               OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), 
                               type=i))
  results$Disease_or_Trait = as.character(results$Disease_or_Trait)
}

results=subset(results, OR>1 & pvalue<0.01 & observed>3)
table(results$type)
#   HTNE    exons promoter   random 
#     82       12       70       2 
results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)

results = results[with(results, order(type, -OR)), ]
write.table(results, paste0("eRNA.SNP.enrichments.",type,".xls"), sep="\t", col.names = T, row.names = F)

results = subset(results, type!='pHTNE') # using private-HTNE

pdf(paste0("eRNA.SNP.enrichments.",type,".pdf"), width=24, height=12); 
# Note: Don't use ggsave() with Rscript, which will generate another Rplot.pdf unnecessarily. See http://stackoverflow.com/questions/19382384/ggplot2-overwrite-one-another-in-rplots-pdf

# re-order the levels in the order of appearance in the data.frame
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=OR, fill=type, ymax=max(OR)*1.1)) + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) + geom_text(aes(label=paste0(observed," (",round(-log10(pvalue),1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by OR)")) 

print(p);

results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by pvalue)")) 

print(p);

dev.off();

# results = results[with(results, order(type, pvalue)), ]
# results = subset(results, type=='HTNE')
# results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
# p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), ymax=max(-log10(pvalue))*1.2)) + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) + ggtitle(paste0("SNP enrichments (LD from ",type,", sorted by pvalue) - HTNE only")) 
# 
# pdf(paste0("eRNA.SNP.enrichments.",type,".pvalue.onlyHTNE.pdf"), width=10, height=12); print(p); dev.off();