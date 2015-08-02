# to test if a TF is more enriched in a region (vs. its background region)
              TF_present TF_absent
DNase_w_eRNA
DNase_wo_eRNA


## all
cut -f1-3 externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 | sortBed | mergeBed -i stdin | awk '{s+=($3-$2);}END{print s}'
# 380,355,257
cut -f4 externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > all.TFBSoccurance.txt
# random
intersectBed -a eRNA.random.bed -b externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 -wo | cut -f4,8 | sort -u | cut -f2 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > random.TFBSoccurance.txt


# eRNA
intersectBed -a eRNA.bed -b externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 -wo | cut -f4,8 | sort -u | cut -f2 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > eRNA.TFBSoccurance.txt

# eRNA w/ DNase peak
intersectBed -b externalData/DNase/merged.DNase.pval.signal.peaks -a eRNA.bed -u | intersectBed -a stdin -b externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 -wo | cut -f4,8 | sort -u | cut -f2 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > eRNAcore.TFBSoccurance.txt

# DNase peak
intersectBed -a externalData/DNase/merged.DNase.pval.signal.peaks -b externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 -wo | cut -f4,14 | sort -u | cut -f2 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > DNASE.TFBSoccurance.txt
# DNase peak w/ eRNA
intersectBed -a externalData/DNase/merged.DNase.pval.signal.peaks -b eRNA.bed -u | intersectBed -a stdin -b externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 -wo | cut -f4,14 | sort -u | cut -f2 | sort | uniq -c | sed -e 's/^[ \t]*//;' | sort -k2,2 > DNASEcore.TFBSoccurance.txt

R

# any TF enriched in eRNAs (vs. the rest of the genome)
setwd("~/eRNAseq")
tf1='all.TFBSoccurance.txt'
tf2='eRNA.TFBSoccurance.txt'
ntf1=51  # eRNA size in million base pair
ntf2=3137  # genome size in million base pair

# all=read.table(tf1); rownames(all)=all[,2]
# x=read.table(tf2); rownames(x)=x[,2]
# df=cbind(x, all[rownames(x),1]); df=df[,-2]; colnames(df)=c('observed','all')
# #df=cbind(df, enrichment.pvalue=apply(df,1,function(x) phyper(x[1]-1,x[2],ntf1-x[2],ntf2,lower.tail=F)))
# results=cbind(TF=rownames(df), df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type="all.eRNA")

tf1='eRNA.TFBSoccurance.txt'
tf2='eRNAcore.TFBSoccurance.txt'
ntf1=5669
ntf2=71469

all=read.table(tf1); rownames(all)=all[,2]
x=read.table(tf2); rownames(x)=x[,2]
df=cbind(x, all[rownames(x),1]); df=df[,-2]; colnames(df)=c('observed','all')
#subset(cbind(df, pvalue=apply(df,1,function(x) phyper(x[1]-1,x[2],ntf2-x[2],ntf1,lower.tail=F))), pvalue<0.01)
results=cbind(TF=rownames(df), df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type="eRNA.eRNA_w_DNase")

tf1='DNASE.TFBSoccurance.txt'
tf2='DNASEcore.TFBSoccurance.txt'
ntf1=6129
ntf2=258989

all=read.table(tf1); rownames(all)=all[,2]
x=read.table(tf2); rownames(x)=x[,2]
df=cbind(x, all[rownames(x),1]); df=df[,-2]; colnames(df)=c('observed','all')
#subset(cbind(df, pvalue=apply(df,1,function(x) phyper(x[1]-1,x[2],ntf2-x[2],ntf1,lower.tail=F))), pvalue<0.01)
results=rbind(results, cbind(TF=rownames(df), df, pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$p.value), OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],x[2]-x[1],ntf1-x[1],ntf2-ntf1-x[2]+x[1]), nrow = 2), alternative='greater')$estimate), type="DNase.DNase_w_eRNA"))

# plot
results=subset(results, pvalue<0.01 & observed>3)
results = results[with(results, order(type, -OR)), ]
table(results$type)
# eRNA.eRNA_w_DNase DNase.DNase_w_eRNA 
# 154                  7 
write.table(results, "eRNA.TF.enrichments.xls", sep="\t", col.names = T, row.names = F)

pdf("~/Dropbox/PDBrainMap/figures/eRNA/eRNA.TF.enrichment.pdf", width = 6,height = 3)
par(mar=c(5,5,2,4), mgp=c(3,.5,.2))
df=subset(results, type=="eRNA.eRNA_w_DNase" & observed>500)
df=df[with(df, order(-observed)), ]
d=barplot(df$observed,  space=.2, names.arg = df$TF, cex.names = 1, cex.axis = 1, las=1, horiz = T, yaxs = "i", xlim=c(0,2200), xlab="Number of DNase-supported HiTNEs")
#text(x=d, y=df$observed, pos=3, offset =0.1, labels=as.character(cut(df$pvalue, breaks=c(0,0.0001,0.001,0.01), include.lowest=T, labels = c('***',"**","*"))), cex=0.6)
#legend("topright",c("* p<0.01","** p<0.001","*** p<0.0001"), bty='n', cex=.8)
dev.off()

par(mar=c(5,5,2,4), mgp=c(3,.5,.2))
df=subset(results, type=="DNase.DNase_w_eRNA")
df=df[with(df, order(-observed)), ]
d=barplot(df$observed, names.arg = df$TF, cex.names = 1, cex.axis = 1, las=1, horiz = T, yaxs = "i", xlim=c(0,2200), xlab="Number of HiTNEs-overlapped DNase peaks")
text(x=d, y=df$observed, pos=3, offset =0.1, labels=as.character(cut(df$pvalue, breaks=c(0,0.0001,0.001,0.01), include.lowest=T, labels = c('***',"**","*"))), cex=0.6)
legend("topright",c("* p<0.01","** p<0.001","*** p<0.0001"), bty='n', cex=.8)
dev.off()
