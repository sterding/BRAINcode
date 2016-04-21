# Rscript for coverage.barplot.pdf (Fig 1c)

setwd("~/neurogen/rnaseq_PD/results/coverage")
# GENCODE meta-exons (v19)
EXON=system("cut -f1-3 $GENOME/Annotation/Genes/exons.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T) # 122000567
CDS=system("cut -f1-3 $GENOME/Annotation/Genes/cds.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T) # 34966072
# intergenic
intergenic=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed | awk '{s+=($3-$2)}END{print s}'",intern = T)  # 650101182
# EXONs
exons=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed -v | intersectBed -a - -b $GENOME/Annotation/Genes/exons.meta.bed | sortBed | mergeBed -i - | awk '{s+=($3-$2)}END{print s}'",intern = T)  # 101113507
# introns
introns=system("intersectBed -a covered.0.05RPM.HCILB_SNDA.bed -b $GENOME/Annotation/Genes/intergenic.bed -v | intersectBed -a - -b $GENOME/Annotation/Genes/exons.meta.bed -v | awk '{s+=($3-$2)}END{print s}'",intern = T) # 1123233053


pdf("coverage.barplot.pdf", paper='us',width=4, height=4)
df=cbind(as.integer(c(EXON,0,0)),as.integer(c(exons,introns,intergenic)))
rownames(df) = c('exons', 'introns', 'intergenic')
colnames(df)=c("GENCODE", "BRAINCODE");
par(mar=c(4,4,2,4))
d=barplot(as.matrix(df), ylim=c(0,3137161264), col=c('#9ecae1','#9ecae1', '#fc9272','#fec44f'), border =NA, axes=F, ylab="Human genome base pairs (in billion)")
text(x=d, y=apply(df,2,sum),pos=3, offset=.2, paste0(round(100*apply(df,2,sum)/3137161264,1),"%"), cex=2)
axis(2, at=c(0:3)*1e9, labels=0:3)
legend("topleft",col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), rownames(df), bty='n', pch=15)

# only show protein-coding exon [JUST FOR CLEMENS]
df=cbind(0,as.integer(c(exons,introns,intergenic)))
df=rbind(as.integer(c(CDS,0)),df)
rownames(df) = c('CDS','exons', 'introns', 'intergenic')
colnames(df)=c("GENCODE", "BRAINCODE");
par(mar=c(4,4,2,4))
d=barplot(as.matrix(df), ylim=c(0,3137161264), col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), border =NA, axes=F, ylab="Human genome base pairs (in billion)")
text(x=d, y=apply(df,2,sum),pos=3, offset=.2, paste0(round(100*apply(df,2,sum)/3137161264,1),"%"), cex=2)
axis(2, at=c(0:3)*1e9, labels=0:3)
legend("topleft",col=c('#3182bd','#9ecae1', '#fc9272','#fec44f'), rownames(df), bty='n', pch=15)

dev.off()