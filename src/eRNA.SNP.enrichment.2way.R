## Rscript to visualize the enrichment of GWAS SNPs in HTNE in a two- or three ways

args<-commandArgs(TRUE)
mergedfile=args[1]

library(tidyverse) # install.packages("tidyverse")  # include ggplot2, tibble, tidyr, readr, purrr, and dplyr. 
library(ggplot2)
library(threejs) # install.packages('threejs')
library(htmlwidgets)
library(reshape2)
library(dplyr)

# setwd("~/eRNAseq"); mergedfile="eRNA.SNP.enrichments.SNAP.merge.wide.qvalue0.01.xls";
CUTOFF=as.numeric(sub(".*qvalue(.*).xls","\\1",mergedfile))

df=read.table(mergedfile, sep="\t", header = T, quote = "\"")

# ===================================
# three-way scatter plot
# ===================================

df %>% mutate_each(funs(ifelse(is.na(.),1,.)), contains("pvalue")) %>% # replace NA to 1
  mutate_each(funs(-log10(.)), contains("pvalue")) %>%
  select(Disease_or_Trait,contains("pvalue")) %>% 
  rename_(.dots=setNames(names(.), gsub("pvalue.", "", names(.)))) %>%
  with(saveWidget(scatterplot3js(HCILB_SNDA,HC_PY,HC_nonNeuron,axisLabels=c("SNDA", "Non-neuron","PY"),
                 grid=T, labels=Disease_or_Trait,flip.y=F,renderer="canvas"),
                 file=sub("xls","html",mergedfile)))

# ===================================
# three-way venn diagram
# ===================================

dd = df %>% select(contains('pvalue')) %>% 
  mutate_each(funs(ifelse(is.na(.),1,.))) %>% 
  mutate_each(funs(.<CUTOFF/1037)) %>%
  mutate(trait=df$Disease_or_Trait) %>% 
  rename_(.dots=setNames(names(.), gsub("pvalue.", "", names(.)))) %>%
  select(trait, HCILB_SNDA, HC_PY, HC_nonNeuron) %>% 
  arrange(desc(HCILB_SNDA), desc(HC_PY), desc(HC_nonNeuron))
  
write.table(dd, file="eRNA.SNP.enrichment.venndiagram.list.xls", quote = T, sep="\t", row.names = F, col.names = T)

dd = dd %>% column_to_rownames(var = "trait") %>% table
d = c(dd)
ns = do.call(paste0, expand.grid(lapply(sub("pvalue..*_","",names(dimnames(dd))),function(x) paste0(x,dimnames(dd)[[1]]))))  # super cool!  # http://stackoverflow.com/a/6254606/951718
names(d) = gsub("TRUE","",gsub("nonNeuronFALSE","",gsub("SNDAFALSE","",gsub("PYFALSE","",ns))))
sink("eRNA.SNP.enrichment.venndiagram.txt")
d
sink()
# run eulerAPE_3 client and make the figure there: eRNA.SNP.enrichment.venndiagram.pdf

# ===================================
# two-way scatter plot
# ===================================

df$Disease_or_Trait = gsub("attention deficit-hyperactivity disorder","ADHD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("Autism spectrum disorder","ASD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("major depressive disorder","MDD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("attention-deficit/hyperactivity disorder","ADHD",df$Disease_or_Trait, ignore.case =T)

pdf(sub("xls","pdf",mergedfile), width=12, height=8); 

# SNDA vs. PY
results = df %>% filter(pmin(pvalue.HCILB_SNDA, pvalue.HC_PY, na.rm =T) <= CUTOFF/1037) %>% 
  reshape(idvar ='Disease_or_Trait',varying=2:ncol(df), sep = ".", # wide->long (http://stackoverflow.com/a/21690303/951718)
          times=c('HCILB_SNDA','HC_nonNeuron','HC_PY'), timevar='type',direction = 'long') %>% 
  na.omit() %>%
  filter(grepl('SNDA|PY',type)) %>%
  mutate(type=factor(as.character(type), levels=c('HCILB_SNDA','HC_PY'))) %>% 
  arrange(type, pvalue) %>%
  mutate(Disease_or_Trait2 = factor(Disease_or_Trait, unique(as.character(Disease_or_Trait)))) 

p=ggplot(results, aes(x=Disease_or_Trait2)) + 
  geom_bar(data = subset(results, type == "HCILB_SNDA"), 
           aes(y=-log10(pvalue), fill = 'SNDA'), width=.5, stat = "identity") +
  geom_text(data = subset(results, type == "HCILB_SNDA"),
            aes(y=-log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2)
p= p + geom_bar(data = subset(results, type == "HC_PY"), 
             aes(y=log10(pvalue), fill = 'PY'), width=.5, stat = 'identity') + 
  geom_text(data = subset(results, type == "HC_PY"),
            aes(y=log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = -90, size=2) + 
  geom_hline(yintercept=c(log10(CUTOFF/1037),-log10(CUTOFF/1037)), size=.5,linetype = 2) +
  scale_x_discrete(drop = FALSE) +  #http://stackoverflow.com/a/15626036/951718 
  ylim(c(log10(min(results$pvalue)),-log10(min(results$pvalue)))*1.1) + 
  xlab("") + ylab("-log10(pvalue)") + theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7), 
        legend.justification=c(1,1), legend.position=c(1,1)) +
  ggtitle("SNDA vs. PY -- SNP enrichments (sorted by pvalue of SNDA)")
print(p);

# SNDA vs. NN
results = df %>% filter(pmin(pvalue.HCILB_SNDA, pvalue.HC_nonNeuron, na.rm =T) <= CUTOFF/1037) %>% 
  reshape(idvar ='Disease_or_Trait',varying=2:ncol(df), sep = ".",   # wide --> long
          times=c('HCILB_SNDA','HC_nonNeuron','HC_PY'), timevar='type',direction = 'long') %>% 
  na.omit() %>%
  filter(grepl('SNDA|non',type)) %>%
  mutate(type=factor(as.character(type), levels=c('HCILB_SNDA','HC_nonNeuron'))) %>% 
  arrange(type, pvalue) %>%
  mutate(Disease_or_Trait2 = factor(Disease_or_Trait, unique(as.character(Disease_or_Trait)))) 

p=ggplot(results, aes(x=Disease_or_Trait2)) + 
  geom_bar(data = subset(results, type == "HCILB_SNDA"), 
           aes(y=-log10(pvalue), fill = 'SNDA'), width=.5, stat = "identity") +
  geom_text(data = subset(results, type == "HCILB_SNDA"),
            aes(y=-log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2)
p= p + geom_bar(data = subset(results, type == "HC_nonNeuron"), 
                aes(y=log10(pvalue), fill = 'nonNeuron'), width=.5, stat = 'identity') + 
  geom_text(data = subset(results, type == "HC_nonNeuron"),
            aes(y=log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = -90, size=2) + 
  geom_hline(yintercept=c(log10(CUTOFF/1037),-log10(CUTOFF/1037)), size=.5,linetype = 2) +
  scale_x_discrete(drop = FALSE) +  #http://stackoverflow.com/a/15626036/951718 
  ylim(c(log10(min(results$pvalue)),-log10(min(results$pvalue)))*1.1) + 
  xlab("") + ylab("-log10(pvalue)") + theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7), 
        legend.justification=c(1,1), legend.position=c(1,1)) +
  ggtitle("SNDA vs. nonNeuron -- SNP enrichments (sorted by pvalue of SNDA)")
print(p);

#  PY vs. NN
results = df %>% filter(pmin(pvalue.HC_PY, pvalue.HC_nonNeuron, na.rm =T) <= CUTOFF/1037) %>% 
  reshape(idvar ='Disease_or_Trait',varying=2:ncol(df), sep = ".", 
          times=c('HCILB_SNDA','HC_nonNeuron','HC_PY'), timevar='type',direction = 'long') %>% 
  na.omit() %>%
  filter(grepl('PY|non',type)) %>%
  mutate(type=factor(as.character(type), levels=c('HC_PY','HC_nonNeuron'))) %>% 
  arrange(type, pvalue) %>%
  mutate(Disease_or_Trait2 = factor(Disease_or_Trait, unique(as.character(Disease_or_Trait)))) 

p=ggplot(results, aes(x=Disease_or_Trait2)) + 
  geom_bar(data = subset(results, type == "HC_PY"), 
           aes(y=-log10(pvalue), fill = 'PY'), width=.5, stat = "identity") +
  geom_text(data = subset(results, type == "HC_PY"),
            aes(y=-log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2)
p= p + geom_bar(data = subset(results, type == "HC_nonNeuron"), 
                aes(y=log10(pvalue), fill = 'nonNeuron'), width=.5, stat = 'identity') + 
  geom_text(data = subset(results, type == "HC_nonNeuron"),
            aes(y=log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = -90, size=2) + 
  geom_hline(yintercept=c(log10(CUTOFF/1037),-log10(CUTOFF/1037)), size=.5,linetype = 2) +
  scale_x_discrete(drop = FALSE) +  #http://stackoverflow.com/a/15626036/951718 
  ylim(c(log10(min(results$pvalue)),-log10(min(results$pvalue)))*1.1) + 
  xlab("") + ylab("-log10(pvalue)") + theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7), 
        legend.justification=c(1,1), legend.position=c(1,1)) +
  ggtitle("PY vs. nonNeuron -- SNP enrichments (sorted by pvalue of PY)")
print(p);

dev.off();

quit('no');







##

# setwd("~/eRNAseq"); mergedfile="eRNA.SNP.enrichments.SNAP.merge.xls";

df=read.table(mergedfile, sep="\t", header = T, quote = "\"")

# ===================================
# three-way scatter plot
# ===================================

df %>% filter(grepl('HTNE',type)) %>% 
  mutate(pvalue=-log10(pvalue)) %>%
  dcast(Disease_or_Trait ~ type, value.var = 'pvalue', fill=0) %>%
  with(saveWidget(scatterplot3js(HTNE.HCILB_SNDA,HTNE.HC_PY,HTNE.HC_nonNeuron,axisLabels=c("SNDA", "Non-neuron","PY"),
                                 grid=T, labels=Disease_or_Trait,flip.y=F,renderer="canvas"),
                  file=sub("xls","html",mergedfile)))

# ===================================
# two-way scatter plot
# ===================================

df$Disease_or_Trait = gsub("attention deficit-hyperactivity disorder","ADHD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("Autism spectrum disorder","ASD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("major depressive disorder","MDD",df$Disease_or_Trait, ignore.case =T)
df$Disease_or_Trait = gsub("attention-deficit/hyperactivity disorder","ADHD",df$Disease_or_Trait, ignore.case =T)

pdf(sub("xls","pdf",mergedfile), width=12, height=8); 

# SNDA vs. PY
results = df %>% filter(grepl('SNDA|PY',type))
results=subset(results, OR>1 & pvalue<0.01/1037 & observed>3)  # only show disease passing the Bonferroni correction
results$type=factor(results$type, levels=c('HTNE.HCILB_SNDA','HTNE.HC_PY'))

results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))

p=ggplot(results, aes(x=Disease_or_Trait2)) + 
  geom_bar(data = subset(results, type == "HTNE.HCILB_SNDA"), 
           aes(y=-log10(pvalue), fill = 'SNDA'), width=.5, stat = "identity") +
  geom_text(data = subset(results, type == "HTNE.HCILB_SNDA"),
            aes(y=-log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5)
p + geom_bar(data = subset(results, type == "HTNE.HC_PY"), 
             aes(y=log10(pvalue), fill = 'PY'), width=.5, stat = 'identity') + 
  geom_text(data = subset(results, type == "HTNE.HC_PY"),
            aes(y=log10(pvalue),label=paste0(observed," (",round(OR,1),")")), 
            position = position_dodge(width=1), hjust=0, vjust=.5, angle = -90, size=2.5) + 
  geom_hline(yintercept=c(log10(0.01/1037),-log10(0.01/1037)), size=.5,linetype = 2) +
  scale_x_discrete(drop = FALSE) +  #http://stackoverflow.com/a/15626036/951718 
  ylim(c(log10(min(results$pvalue)),-log10(min(results$pvalue)))*1.1) + 
  xlab("") + ylab("-log10(pvalue)") + theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        legend.justification=c(1,1), legend.position=c(1,1)) +
  ggtitle("SNDA vs. PY -- SNP enrichments (sorted by pvalue of SNDA)")


p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle("SNDA vs. PY -- SNP enrichments (sorted by pvalue)")
#p = p + coord_flip() 
print(p);

# SNDA vs. NN
results = df %>% filter(grepl('SNDA|non',type))
results=subset(results, OR>1 & pvalue<0.01/1037 & observed>3)  # only show disease passing the Bonferroni correction
results$type=factor(results$type, levels=c('HTNE.HCILB_SNDA','HTNE.HC_nonNeuron'))

results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle("SNDA vs. NN -- SNP enrichments (sorted by pvalue)")
#p = p + coord_flip() 
print(p);

#  PY vs. NN
results = df %>% filter(grepl('PY|non',type))
results=subset(results, OR>1 & pvalue<0.01/1037 & observed>3)  # only show disease passing the Bonferroni correction
results$type=factor(results$type, levels=c('HTNE.HC_PY','HTNE.HC_nonNeuron'))

results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle("PY vs. NN -- SNP enrichments (sorted by pvalue)") 
#p = p + coord_flip() 
print(p);

dev.off();

