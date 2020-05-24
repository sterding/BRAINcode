###########################################
# Rscript to plot the expression of specific gene(s) along the stages/categories
# Usage: Rscript $PIPELINE_PATH/_plotTrend.R SNCA.fpkm.allSamples.uniq.xls SNCA.pdf
# or
# grep -w -P "tracking_id|SNCA|GBA|LRRK2" genes.fpkm.allSamples.uniq.xls | Rscript ~/neurogen/pipeline/RNAseq/modules/_plotTrend.R stdin SCNA.pdf
# grep -w -P "tracking_id|SNCA|GBA|LRRK2" SRR_tissue.genes.fpkm.allSamples.uniq.xls | Rscript ~/neurogen/pipeline/RNAseq/modules/_plotTrend.R stdin SNCA.pdf GTEx
# Author: Xianjun Dong
# Version: 0.0
# Date: 2014-Jul-15
###########################################

args<-commandArgs(TRUE)

require(ggplot2)
require(grid)
source("~/neurogen/pipeline/RNAseq/bin/geomBoxplot_noOutliers.R")


FPKMfile=args[1]  # either filename or stdin
outputfile=args[2]
projectCODE=ifelse(is.na(args[3]),"BRAINCODE",args[3])  # GTEx or PD

fpkm=read.table(file(FPKMfile), sep="\t", header=T);  # table with header (1st row) and ID (1st column)
# change from wide to long format
library('reshape2')

if(projectCODE=="GTEx"){
    df=melt(fpkm, measure.vars=grep("SRR", colnames(fpkm)))
    df$variable=gsub("SRR[0-9]*.(.*)", "\\1", df$variable) 
} else{
    colnames(fpkm) = gsub("FPKM.","",colnames(fpkm))
    df=melt(fpkm, measure.vars=grep("_rep", colnames(fpkm)))
    df$variable=gsub("(.*)_.*_(.*)_.*_.*", "\\1_\\2", df$variable)
}

df$tracking_id=paste(df$tracking_id, " (", df$gene_short_name, ")", sep="")

# plot
pdf(file=outputfile, width = 10, height = 5)
p=ggplot(df, aes(colour=tracking_id, y=value, x=variable))+
  geom_boxplot_noOutliers() +
  #geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0,50)) + 
  ylab("FPKM") +
  theme(legend.justification=c(0,1),
        legend.position=c(0,1),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size =  unit(.4, "cm"),
        axis.text.x=element_text(angle=-90, hjust = 0)
  )
print(p)


p = ggplot(df, aes(y=value, x=variable)) + 
  geom_boxplot() + 
  theme_bw() + 
  facet_wrap(~tracking_id, scales = "free")

print(p)
dev.off()

##fpkm=read.table("genes.fpkm.allSamples.uniq.xls", header=T)
#fpkm=cbind(id=paste(fpkm$tracking_id, fpkm$gene_id, fpkm$gene_short_name, sep="__"), fpkm[, grep("FPKM", colnames(fpkm))])
##rownames(fpkm)=fpkm[,1]; fpkm=fpkm[,-1]
#
## change from wide to long format
#library('reshape2')
#fpkml=melt(fpkm)
#df=fpkml[grep(pattern,fpkml$id), ]
#df$variable=gsub("FPKM.(.*)_.*_(.*)_.*", "\\1_\\2", df$variable)
#
## plot
#library('ggplot2')
#p <- ggplot(df, aes(colour=id, y=value, x=variable))
#p +  geom_boxplot() +
#     theme(legend.justification=c(0,1),
#           legend.position=c(0,1),
#           legend.title=element_blank(),
#           legend.background = element_rect(fill="transparent")
#        )
#
##p + geom_point() + geom_errorbar(aes(ymax = resp + se, ymin=resp - se), width=0.2)
##p + geom_line(aes(group=group)) + geom_errorbar(limits, width=0.2)


