## Rscript to draw the venPier plot
# see http://onetipperday.blogspot.com/2014/09/vennpier-combination-of-venn-diagram.html
# Usage:
# Rscript $pipeline_path/modules/_bam2annotation.r accepted_hits.bam.bam2annotation accepted_hits.bam.bam2annotation.pdf

args<-commandArgs(TRUE)

stat_file=args[1]  # for example: stat_file="/PHShome/xd010/neurogen/rnaseq_PD/run_output/PD_BN13-18_SNDA_5_rep1/uniq/accepted_hits.bam.bam2annotation"
pdf_file=args[2]

df=read.table(stat_file)
d=as.list(df[,2]); names(d)=gsub(":.*","", gsub("-","_",df[,1])); attach(d);

## internal variables: exons, intergenic, intergenic_notneargene, introns, LINE, mtRNA, rRNA, SINE, total, utr3, utr5, total_non_rRNA_mt
intergenic_not_near_genes=intergenic_notneargene

rest=total_non_rRNA_mt
genic=rest-intergenic
introns_and_exons=introns+exons-genic
intergenic_near_genes = intergenic - intergenic_not_near_genes
 
# parameter for pie chart
iniR=0.2 # initial radius
colors=list(NO='white',total='black', mtRNA='#e5f5e0',rRNA='#a1d99b',genic='#3182bd',exons='#9ecae1',introns='#fc9272',intergenic='#fec44f',intergenic_near_genes='#fee0d2',intergenic_not_near_genes='#d95f0e')
 
#if(!(require('plotrix'))) install.packages('plotrix')  # not working for Rscript
library('plotrix')
 
pdf(pdf_file, width=12, height=6)
par(mfrow=c(1,2))
# from outer circle to inner circle
#0 circle: blank
pie(1, radius=iniR, init.angle=90, col=c('white'), border = NA, labels='')
 
#4 circle: show genic:exons and intergenic:downstream
floating.pie(0,0,c(exons, total-exons),radius=5*iniR, startpos=pi/2, col=as.character(colors[c('exons','NO')]),border=NA)
 
#3 circle: show genic:introns and intergenic:intergenic_notneargene | upstream
floating.pie(0,0,c(genic-introns, introns, intergenic_not_near_genes, intergenic_near_genes, mtRNA+rRNA),radius=4*iniR, startpos=pi/2, col=as.character(colors[c('NO','introns','intergenic_not_near_genes','intergenic_near_genes','NO')]),border=NA)
 
#2 circle: divide the rest into genic and intergenic
floating.pie(0,0,c(genic, intergenic, rRNA, total-rest-rRNA),radius=3*iniR, startpos=pi/2, col=as.character(colors[c('genic','intergenic','rRNA','NO')]),border=NA)
 
#1 circle: for rRNA+mtRNA+rest
floating.pie(0,0, c(total-mtRNA,mtRNA), radius=2*iniR, startpos=pi/2, col=as.character(colors[c('NO','mtRNA')]), border = NA)
 
#legend(0, 5*iniR, gsub("intergenic ", "--", gsub("_"," ",names(colors)[-1])), col=as.character(colors[-1]), pch=19,bty='n', ncol=2, cex=0.6)

names=gsub("intergenic ", "--", gsub("_"," ",names(colors)[-1]))
values = sapply(names(colors)[-1], get)
percent=format(100*values/total, digits=2, trim=T)
values = format(values, big.mark=",", scientific=FALSE, trim=T)
cl=as.character(colors[-1])
pchs=rep(19, length(cl)); pchs[1]=1;
legend(0, 5*iniR, paste(names," (",values,", ", percent,"%)", sep=""), col=cl, pch=pchs,bty='n', ncol=1, cex=0.6)


## a sub vennpieR plot only for the rest reads
# =========================================
pie(1, radius=iniR, init.angle=90, col=c('white'), border = NA, labels='')
floating.pie(0,0,c(exons, rest-exons),radius=3*iniR, startpos=pi/2, col=as.character(colors[c('exons','NO')]),border=NA)
floating.pie(0,0,c(genic-introns, introns, intergenic_not_near_genes, intergenic_near_genes),radius=2*iniR, startpos=pi/2, col=as.character(colors[c('NO','introns','intergenic_not_near_genes','intergenic_near_genes')]),border=NA)
floating.pie(0,0,c(genic, intergenic),radius=iniR, startpos=pi/2, col=as.character(colors[c('genic','intergenic')]),border=NA)

colors=c(rest="black", colors[-(1:4)])
names=gsub("intergenic ", "--", gsub("_"," ",names(colors)))
values = sapply(names(colors), get)
percent=format(100*values/rest, digits=2, trim=T)
values = format(values, big.mark=",", scientific=FALSE, trim=T)
cl=as.character(colors)
pchs=rep(19, length(cl)); pchs[1]=1;
legend(0, 5*iniR, paste(names," (",values,", ", percent,"%)", sep=""), col=cl, pch=pchs,bty='n', ncol=1, cex=0.6)

dev.off()