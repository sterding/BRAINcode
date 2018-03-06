## Rscript to test the enrichment of GWAS SNPs in HTNE
#
# The Fisher exact test was performed to test the odds ratio of trait association with enhancers vs non-enhancers taking into account the genomic distribution of SNPs.  
# OR = (A/B) / (C/D), where
# 
# A: number of trait-associated SNPs in enhancers
# B: number of non-trait-associated SNPs in enhancers (from dbSNP - GWAS SNPs for that trait)
# C: number of trait-associated SNPs distal to enhancers
# D: number of non-trait-associated SNPs distal to enhancers (from dbSNP - GWAS SNPs for that trait)

args<-commandArgs(TRUE)
type=args[1]
samplegroup=args[2]
output_prefix=args[3]

if(length(args)==2) output_prefix=paste0("eRNA.SNP.enrichments.",type)

require(tidyverse)

setwd(paste0("~/eRNAseq/",samplegroup));

# setwd("~/eRNAseq/HCILB_SNDA"); type="PICS"; output_prefix=paste0("eRNA.SNP.enrichments.",type)

df=read.table(paste0("SNP.",type,".counts.summary"), header=F, col.names = c('type',"N_SNP","BP_length"), stringsAsFactors=F); 
df = df %>% mutate(SNP_density = N_SNP*1e6/BP_length)
selected_to_display=c('total','coding_exon','intron','intergenic','promoter','FANTOM5_enhancer','chromHMM_enhancer','chromHMM_enhancer_brain','TNE','random')
df = df %>% filter(type %in% selected_to_display)
df$type = factor(df$type, rev(selected_to_display))

pdf(paste0(output_prefix,".barplot.pdf"), width = 5, height = 4)
p = ggplot(df, aes(x=type, y=SNP_density)) 
p = p + geom_col(width = 0.8,fill='#4EB3CD') +  coord_flip()
p = p + theme_bw() 
#p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=N_SNP), col='white',hjust=1, vjust=.5) 
p = p + labs(title=paste0(basename(getwd()), " -- ",type," SNPs"), x="", y="SNP density (SNPs/million base pair)") 
print(p);
dev.off()