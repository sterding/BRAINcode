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

require(ggplot2)

setwd(paste0("~/eRNAseq/",samplegroup));

# setwd("~/eRNAseq/HCILB_SNDA"); type="SNAP"; output_prefix=paste0("eRNA.SNP.enrichments.",type)
# setwd("~/eRNAseq/HC_PY"); type="SNAP"; output_prefix=paste0("eRNA.SNP.enrichments.",type)
# setwd("~/eRNAseq/HC_nonNeuron"); type="SNAP"; output_prefix=paste0("eRNA.SNP.enrichments.",type)

s=read.table(paste0("SNP.",type,".counts.summary"), header=F,row.names=1); 
results=data.frame();
for(i in c('HTNE','HTNE1','HTNE2','HTNE3','HTNE-private','exon','promoter','random','chromHMM_cellline','chromHMM_brain','DNase','CAGE','TFhotspot','P300','HCNE')){
  if(!file.exists(paste0("SNP.",type,".count.",i))) next;
  n1=s[i,1]; n2=s['all',1];
  all=read.table(paste0("SNP.",type,".count.all")); rownames(all)=all[,1]
  x=read.table(paste0("SNP.",type,".count.",i)); rownames(x)=x[,1]
  df=cbind(x, all[rownames(x),2]); # only the traits occurred
  #df=merge(all, x, by="V1", all=T); df[is.na(df)]=0; rownames(df)=df[,1]; # all traits (0 for nonoccurance)  --> not working for fisher.test
  df=df[,-1]; colnames(df)=c('observed','all')
  results=rbind(results, cbind(Disease_or_Trait=rownames(df), 
                               df, 
                               pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), 
                               OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1], n1-x[1], n2-n1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), 
                               type=i))
  results$Disease_or_Trait = as.character(results$Disease_or_Trait)
}

results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
results = results[with(results, order(type, -OR)), ]

# save all result
write.table(results, paste0(output_prefix,".full.xls"), sep="\t", col.names = T, row.names = F)

results=subset(results, OR>1 & pvalue<0.01/1037 & observed>3)
table(results$type)
# HTNE             HTNE1             HTNE2             HTNE3              exon          promoter            random chromHMM_cellline    chromHMM_brain 
# 61                26                20                49                11                43                 0               264               151 
# DNase              CAGE         TFhotspot              P300              HCNE 
# 115                36               201                65                 0 

# add customerized category info (defined by Clemens, annotated till HTNE, HTNE1, HTNE2, HTNE3, exon, promoter, random, chromHMM) 
traits_dopamine="response to antipsychotic treatment|sleep quality|insomnia|sleep-related phenotypes|schizophrenia|Parkinson|attention deficit hyperactivity disorder|response to methylphenidate treatment|addiction|bipolar disorder|^Response to iloperidone treatment|dyskinesia|post-traumatic stress disorder|Nicotine dependence|Major mood disorders|Narcolepsy|Narcolepsy with cataplexy|Tourette syndrome|Psychosis and Alzheimer|Alcohol consumption|Restless legs syndrome"
traits_heart="cardiac|QRS duration|heart|coronary|echocardiographic|tetralogy of fallot|PR interval|QT interval|Abdominal aortic aneurysm|P wave duration"
traits_obesity="underweight|insulin resistance|waist circumference|diabetes|IgG glycosylation|obesity|cholesterol|Metabolic|Very long-chain saturated fatty acid levels|Proinsulin levels|Red blood cell fatty acid levels|Trans fatty acid levels|Glycated hemoglobin levels|Adiponectin levels|Body mass index in asthmatics|Insulin-like growth factors|Height adjusted BMI|Glycemic traits|Triglycerides-Blood Pressure|Response to fenofibrate|Response to statin therapy|Plasma omega-6 polyunsaturated fatty acid levels|Phospholipid levels|Waist-hip ratio|Triglycerides|Fasting glucose-related traits|IgG glycosylation|Body mass index|Obesity|Comprehensive strength and appendicular lean mass|Lipoprotein|Body mass in chronic obstructive pulmonary disease"
traits_brainvolume="brain region volumes|MRI atrophy measure|grey matter volume|Normalized brain volume"
results$Trait_category = ifelse(grepl(traits_dopamine,results$Disease_or_Trait, ignore.case = T),"Dopamine system traits",ifelse(grepl(traits_heart,results$Disease_or_Trait, ignore.case = T),"Cardiovascular traits",ifelse(grepl(traits_obesity,results$Disease_or_Trait, ignore.case = T),"Obesity, weight, diabetes traits", ifelse(grepl(traits_brainvolume,results$Disease_or_Trait, ignore.case = T),"Brain volume traits","Other"))))

results$Trait_category=factor(results$Trait_category, levels=rev(c('Dopamine system traits','Cardiovascular traits','Obesity, weight, diabetes traits',"Brain volume traits","Other")))

write.table(results, paste0(output_prefix,".qvalue0.01.annotated.xls"), sep="\t", col.names = T, row.names = F)

## plot barplot of total traits per category
## ==========================================

# setwd("~/eRNAseq/HCILB_SNDA"); type="SNAP"; output_prefix=paste0("eRNA.SNP.enrichments.",type)
# results=read.table(paste0(output_prefix,".qvalue0.01.annotated.xls"), header = T, sep = "\t"); head(results)
# results$type = factor(results$type, levels=c('CAGE', 'chromHMM_cellline','chromHMM_brain', 'DNase', 'exon', 'HTNE', 'HTNE1', 'HTNE2', 'HTNE3', 'P300', 'promoter', 'TFhotspot','random'))

pdf(paste0(output_prefix,".barplot.pdf"), paper = 'us')
par(mar=c(4,6,4,2))

## TNE/exon/promoter/random
barplot(t(table(results$type, results$Trait_category)[c('HTNE','exon','promoter','random'),]), col=rev(c('deeppink3','deepskyblue','orange','forestgreen','gray')), names.arg=c('TNE\nsite','Exon','Promoter','Random\nbackground'), ylab="Number of GWAS traits/diseases\nwith over-represented disease variants", cex.names=0.7, space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))

## TNE and subclass
df=table(results$type, results$Trait_category)[c('HTNE','HTNE1', 'HTNE2', 'HTNE3'),]
barplot(t(df), col=rev(c('deeppink3','deepskyblue','orange','forestgreen','gray')), names.arg=c('TNE','TNE (class I)', 'TNE (class II)', 'TNE (class III)'), ylab="Number of GWAS traits/diseases\nwith over-represented disease variants", cex.names=0.7, space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))
# density
xx=barplot(t(100*df/rowSums(df+0.0001)), col=rev(c('deeppink3','deepskyblue','orange','forestgreen','gray')), names.arg=c('TNE','TNE (class I)', 'TNE (class II)', 'TNE (class III)'), ylab="Density of GWAS traits/diseases\nwith over-represented disease variants (%)", cex.names=0.7, ylim=c(0,100),space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))
mtext(text = round(100*df/rowSums(df+0.0001),1)[,'Dopamine system traits'], side=1, line=-1.2, at=xx, cex = .8)

## TNE/CAGE/chromHMM/DNase/P300 and denstiy
df=table(results$type, results$Trait_category)[c('HTNE','CAGE','chromHMM_cellline','DNase','P300'),]
# density only for DA
xx=barplot(t(100*df[,'Dopamine system traits']/rowSums(df+0.0001)), col='gray', names.arg=c('Dopamine neuron\nTNE sites','CAGE\nenhancers','chromHMM\nenhancers','DHS\nsites','P300\npeaks'), ylab="Percentage of Dopamine system traits\nwith over-represented disease variants (%)", ylim=c(0,20),cex.names=0.7,space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))
mtext(text = df[,'Dopamine system traits'], side=1, line=-1.2, at=xx, cex = .8)
# density
xx=barplot(t(100*df/rowSums(df+0.0001)), col=rev(c('deeppink3','deepskyblue','orange','forestgreen','gray')), names.arg=c('Dopamine neuron\nTNE sites','CAGE\nenhancers','chromHMM\nenhancers','DHS\nsites','P300\npeaks'), ylab="Density of GWAS traits/diseases\nwith over-represented disease variants (%)", cex.names=0.7, ylim=c(0,100),space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))
mtext(text = rowSums(df), side=1, line=-1.2, at=xx, cex = .8)
# raw count
barplot(t(df), col=rev(c('deeppink3','deepskyblue','orange','forestgreen','gray')), names.arg=c('Dopamine neuron\nTNE sites','CAGE\nenhancers','chromHMM\nenhancers','DHS\nsites','P300\npeaks'), ylab="Number of GWAS traits/diseases\nwith over-represented disease variants", cex.names=0.7, space=0.9,border = NA, axis.lty=1, legend.text=T, args.legend = list(x='topright',bty = "n", cex=.7, border =NA))

dev.off()


pdf(paste0(output_prefix,".pdf"), width=20, height=12); 
# Note: Don't use ggsave() with Rscript, which will generate another Rplot.pdf unnecessarily. See http://stackoverflow.com/questions/19382384/ggplot2-overwrite-one-another-in-rplots-pdf

# re-order the levels in the order of appearance in the data.frame

results=subset(results, pvalue<0.01/1037)  # only show disease passing the Bonferroni correction

results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by pvalue)")) 

print(p);

results = results[with(results, order(type, -OR)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=OR, fill=type, ymax=max(OR)*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(-log10(pvalue),1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by OR)")) 

print(p);

dev.off();

## only HTNE, HTNE1, 2, and 3
pdf(paste0(output_prefix,".HTNE123.pdf"), width=20, height=12); 
results = subset(results, type %in% c('HTNE','HTNE1','HTNE2','HTNE3'))
results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by pvalue)")) 

print(p);

dev.off();

## private HTNE only
pdf(paste0(output_prefix,".HTNE-private.pdf"), width=20, height=12); 
results = subset(results, type=='HTNE-private')
results = results[with(results, order(type, pvalue)), ]
results$Disease_or_Trait2 <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
p = ggplot(results, aes(x=Disease_or_Trait2, y=-log10(pvalue), fill=type, ymax=max(-log10(pvalue))*1.1)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + geom_hline(yintercept=-log10(0.01/1037), size=.5,linetype = 2)  ## Bonferroni correction, where 1037 is the number of disease/traits in GWAS (wc -l SNP.$type.count.all)
p = p + theme_bw() 
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.justification=c(1,1), legend.position=c(1,1)) 
p = p + geom_text(aes(label=paste0(observed," (",round(OR,1),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle(paste0(basename(getwd()), " -- SNP enrichments (LD from ",type,", sorted by pvalue)")) 

print(p);

dev.off();

