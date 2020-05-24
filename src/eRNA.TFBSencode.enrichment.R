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
input=args[1]

require(tidyverse)
require(grid)
library(parallel)

# setwd("~/eRNAseq/HC_nonNeuron"); input="eRNA.JASPARmotifscan.2018.txt"

df=read.table(input, header=F); 
colnames(df) = c('ID','AB','AnB', 'nAB', 'nAnB','type');

# results=cbind(df,
# 				pvalue = apply(df[,2:5], 1, function(x) fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow = 2), alternative='greater')$p.value), 
# 				OR=apply(df[,2:5], 1, function(x) fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow = 2), alternative='greater')$estimate)
# 			)

results=cbind(df, do.call("rbind", 
                          apply(df[,2:5], 1, function(x) fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow = 2), alternative='greater')[c('p.value','estimate')])
                          )
              )
names(results)[names(results) == 'p.value'] <- 'pvalue'
names(results)[names(results) == 'estimate'] <- 'OR'
results$pvalue=as.numeric(results$pvalue)
results$OR=as.numeric(results$OR)

N = length(unique(as.character(results$ID)))
results=subset(results, OR>1 & pvalue<0.01/N & AB>3)

# add the minimal non-zero pvalue to all pvalues, to avoid log10(0)
results$pvalue = results$pvalue + min(results$pvalue[results$pvalue>0])

results = results[with(results, order(type, -OR)), ]
write.table(results, paste0(input,".qvalue0.01.xls"), sep="\t", col.names = T, row.names = F)

# setwd("~/eRNAseq/HCILB_SNDA"); input="eRNA.JASPARmotifscan.2018.txt"
# results=read.table(paste0(input,".qvalue0.01.xls"), header = T, sep = "\t", stringsAsFactors = F); N=579; head(results);

pdf(paste0(input,".pdf"), width=7.2, height=3); 
# Note: Don't use ggsave() with Rscript, which will generate another Rplot.pdf unnecessarily. See http://stackoverflow.com/questions/19382384/ggplot2-overwrite-one-another-in-rplots-pdf

#results$pvalue[results$pvalue<1e-50]=1e-50
#results$OR[results$OR>20]=20

result = results[with(results, order(type, pvalue)), ]
#result = head(result, 50) # top 50 only
result$ID <- factor(result$ID, unique(as.character(result$ID)))
p = ggplot(result, aes(x=ID, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/N), size=.5,linetype = 2)  ## Bonferroni correction, where N is the number of TF motifs in JASPAR
p = p + theme_bw() 
p = p + theme(plot.title = element_text(size = 5), axis.title.y = element_text(size=5),axis.text.y = element_text(size=5), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=5), legend.text=element_text(size = 5), legend.key.size=unit(5,'pt'),legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(AB," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=1.5) 
p = p + ggtitle(paste0(basename(input),", sorted by pvalue)")) 

print(p);

result = results[with(results, order(type, -OR)), ]
#result = head(result, 50) # top 50 only
result$ID <- factor(result$ID, unique(as.character(result$ID)))
p = ggplot(results, aes(x=ID, y=OR, fill=type, ymax=max(OR)*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + theme_bw() 
p = p + theme(plot.title = element_text(size = 5), axis.title.y = element_text(size=5),axis.text.y = element_text(size=5), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=5), legend.text=element_text(size = 5), legend.key.size=unit(5,'pt'), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(AB," (",round(-log10(pvalue),1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=1.5) 
p = p + ggtitle(paste0(basename(input),", sorted by OR)")) 

print(p);

## each type separately
if(length(unique(as.character(results$type)))>1){
  for(tp in unique(as.character(results$type))){
    result = subset(results, type==tp) %>% arrange(pvalue) %>% head(50)  # top 50 only
    result$ID <- factor(result$ID, unique(as.character(result$ID)))
    result$type <- as.character(result$type)
    p = ggplot(result, aes(x=ID, y=-log10(pvalue), ymax=max(-log10(pvalue))*1.1)) 
    p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
    p = p + geom_hline(yintercept=-log10(0.01/N), size=.5,linetype = 2)  ## Bonferroni correction, where N is the number of TF motifs in JASPAR
    p = p + theme_bw() 
    p = p + theme(plot.title = element_text(size = 5), axis.title.y = element_text(size=5),axis.text.y = element_text(size=5), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=5), legend.text=element_text(size = 5), legend.key.size=unit(5,'pt'),legend.justification=c(1,1), legend.position=c(1,1)) 
    p = p + geom_text(aes(label=paste0(AB," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=1.5) 
    p = p + ggtitle(paste0(basename(input),"- sorted by pvalue - ", tp)) 
    
    print(p);
  }
}

dev.off()