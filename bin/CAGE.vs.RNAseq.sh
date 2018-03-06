# CAGE vs. RNAseq
#!/bin/bash

cd ~/neurogen/CAGE_PDBrainMap/results

CAGEfwd=~/neurogen/CAGE_PDBrainMap/output_dir/HC_UWA616_SN_1_rep1/accepted_hits.plus.normalized.bw
CAGErev=~/neurogen/CAGE_PDBrainMap/output_dir/HC_UWA616_SN_1_rep1/accepted_hits.minus.normalized.bw
RNAseq=~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.accepted_hits.normalized.bw

# get mean RPM for promoter regions of protein-coding Tx
awk '{OFS="\t"; if($4~/protein_coding\.protein_coding/) {split($4,a,"___"); print a[3];}}' $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 > tmp
fgrep -wf tmp ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.isoforms.fpkm_tracking | awk '$13=="OK"' | cut -f1 > pcTx.list # N=67638 

fgrep -f pcTx.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($6=="+") {tss=$2; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.promoter.tab
fgrep -f pcTx.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($6=="-") {tss=$2; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.promoter.tab
fgrep -wf pcTx.list ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.isoforms.fpkm_tracking | awk '$13=="OK"' | cut -f1,10 > RNAseq.SNDA.Tx.tab

# get mean RPM for permissive enhancers
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.cageEnhancer.tab
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.cageEnhancer.tab
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $RNAseq  stdin /dev/null -minMax -bedOut=RNAseq.SNDA.cageEnhancer.tab

# get mean RPM for HiTNE regions
bigWigAverageOverBed $CAGEfwd ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=CAGE.SNDA.fwd.HiTNE.tab
bigWigAverageOverBed $CAGErev ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=CAGE.SNDA.rev.HiTNE.tab
bigWigAverageOverBed $RNAseq  ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=RNAseq.SNDA.HiTNE.tab

# get mean RPM for eHiTNE regions
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SNDA.eHiTNE.tab


CAGEfwd=~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.SN.fwd.bigwig
CAGErev=~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.SN.rev.bigwig
RNAseq=~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.accepted_hits.normalized.bw

fgrep -wf tmp ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.isoforms.fpkm_tracking | awk '$13=="OK"' | cut -f1 > pcTx.list

fgrep -f pcTx.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($6=="+") {tss=$2; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.promoter.tab
fgrep -f pcTx.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($6=="-") {tss=$2; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.promoter.tab
fgrep -wf pcTx.list ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.isoforms.fpkm_tracking | awk '$13=="OK"' | cut -f1,10 > RNAseq.SN.Tx.tab

# get mean RPM for permissive enhancers
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.cageEnhancer.tab
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.cageEnhancer.tab
grep -v track ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $RNAseq  stdin /dev/null -minMax -bedOut=RNAseq.SN.cageEnhancer.tab

# get mean RPM for HiTNE regions
bigWigAverageOverBed $CAGEfwd ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=CAGE.SN.fwd.HiTNE.tab
bigWigAverageOverBed $CAGErev ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=CAGE.SN.rev.HiTNE.tab
bigWigAverageOverBed $RNAseq  ~/eRNAseq/HCILB_SNDA/eRNA.bed /dev/null -minMax -bedOut=RNAseq.SN.HiTNE.tab

# get mean RPM for eHiTNE regions
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($28<3)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/HCILB_SNDA/eRNA.characterize.xls | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SN.eHiTNE.tab

## merge together

join -1 1 -2 1 <(cat CAGE.SNDA.fwd.promoter.tab CAGE.SNDA.rev.promoter.tab | cut -f4,7 | sort) <(sort RNAseq.SNDA.Tx.tab) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","promoter"}'> CAGE.RNAseq.tab
#join -1 1 -2 1 <(cat CAGE.SNDA.fwd.promoter.tab CAGE.SNDA.rev.promoter.tab | cut -f4,7 | sort) <(cut -f4,7 RNAseq.SNDA.promoter.tab | sort) | awk '{OFS="\t"; print $1,($2>0)?$2:(-$2),$3,"SNDA","promoter"}'> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.cageEnhancer.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.cageEnhancer.tab | sort) | awk '{print $1,$2+$3;}') <(cut -f4,7 RNAseq.SNDA.cageEnhancer.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","cageEnhancer"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.HiTNE.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.HiTNE.tab | sort) | awk '{print $1,$2+$3;}') <(cut -f4,7 RNAseq.SNDA.HiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","HiTNE"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.eHiTNE.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.eHiTNE.tab | sort) | awk '{print $1,$2+$3;}') <(cut -f4,7 RNAseq.SNDA.eHiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","eHiTNE"}' >> CAGE.RNAseq.tab

#join -1 1 -2 1 <(cat CAGE.SN.fwd.promoter.tab CAGE.SN.rev.promoter.tab | cut -f4,7 | sort) <(cut -f4,7 RNAseq.SN.promoter.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","promoter"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(cat CAGE.SN.fwd.promoter.tab CAGE.SN.rev.promoter.tab | cut -f4,7 | sort) <(sort RNAseq.SN.Tx.tab) | awk '{OFS="\t"; print $1,$2,$3,"SN","promoter"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.cageEnhancer.tab | sort) <(cut -f4,7 CAGE.SN.rev.cageEnhancer.tab | sort) | awk '{print $1,$2+$3}') <(cut -f4,7 RNAseq.SN.cageEnhancer.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","cageEnhancer"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.HiTNE.tab | sort) <(cut -f4,7 CAGE.SN.rev.HiTNE.tab | sort) | awk '{print $1,$2+$3}') <(cut -f4,7 RNAseq.SN.HiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","HiTNE"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.eHiTNE.tab | sort) <(cut -f4,7 CAGE.SN.rev.eHiTNE.tab | sort) | awk '{print $1,$2+$3}') <(cut -f4,7 RNAseq.SN.eHiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","eHiTNE"}' >> CAGE.RNAseq.tab

##
R
df=read.table("CAGE.RNAseq.tab", header=F,stringsAsFactors=F)
colnames(df)=c("id","CAGE","RNAseq","tissue_regio_for_CAGE","annotation_regions")
#dim(subset(df, CAGE==0 & RNAseq==0))
#df=subset(df, CAGE>0 & RNAseq>0)
pdf("CAGE.RNAseq.pdf", width=10, height=3)
par(mfrow=c(1,4), pty="s"); # to make sure the frame is a square (width=height)

# CAGE (SN) vs. RNAseq (SN)
j='promoter'
ID=intersect(df$id[df$tissue_regio_for_CAGE=='SNDA' & df$annotation_regions==j], df$id[df$tissue_regio_for_CAGE=='SN' & df$annotation_regions==j])
plot(0.01+subset(df, id %in% ID & tissue_regio_for_CAGE=='SN' & annotation_regions==j, select = RNAseq)[,1], 0.01+subset(df, id %in% ID & tissue_regio_for_CAGE=='SNDA' & annotation_regions==j, select = CAGE)[,1], xlab="RNAseq (SN)", ylab="CAGE (SN)", log='xy', main='HC_UWA616', pch=19, cex=0.6, col='#00000044')

# CAGE (SN) vs. RNAseq (SNDA)
i='SNDA'; j='promoter';
plot(CAGE~RNAseq, data=0.01+subset(df, tissue_regio_for_CAGE==i & annotation_regions==j, select=c(CAGE, RNAseq)), log='xy', xlab="RNAseq (SNDA)", ylab="CAGE (SN)", main='HC_UWA616', pch=19, cex=0.6, col='#00000044')

# braincode CAGE vs. fantom CAGE (SN)
j='promoter'
ID=intersect(df$id[df$tissue_regio_for_CAGE=='SN' & df$annotation_regions==j], df$id[df$tissue_regio_for_CAGE=='SNDA' & df$annotation_regions==j])
plot(0.01+subset(df, id %in% ID & tissue_regio_for_CAGE=='SN' & annotation_regions==j, select = CAGE)[,1], 0.01+subset(df, id %in% ID & tissue_regio_for_CAGE=='SNDA' & annotation_regions==j, select = CAGE)[,1], xlab="FANTOM CAGE (SN)", ylab="BRAINCODE CAGE (SN)", log='xy', main='HC_UWA616', pch=19, cex=0.6, col='#00000044')

# RNAseq SN vs. SNDA
j='promoter'
ID=intersect(df$id[df$tissue_regio_for_CAGE=='SN' & df$annotation_regions==j], df$id[df$tissue_regio_for_CAGE=='SNDA' & df$annotation_regions==j])
plot(0.001+subset(df, id %in% ID & tissue_regio_for_CAGE=='SN' & annotation_regions==j, select = RNAseq)[,1], 0.001+subset(df, id %in% ID & tissue_regio_for_CAGE=='SNDA' & annotation_regions==j, select = RNAseq)[,1], xlab="RNAseq (SN)", ylab="RNAseq (SNDA)", log='xy', main='HC_UWA616', pch=19, cex=0.6, col='#00000044')

dev.off()

## SN and SNDA private genes
