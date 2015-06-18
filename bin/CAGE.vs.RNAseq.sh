# CAGE vs. RNAseq
#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Script to generate bigwig files for CAGE bam input"
  echo "=================================================="
  echo "Usage: $0 input.bam"
  echo "Output: input.plus.bw and input.minus.bw"
  exit
fi

CAGEfwd=~/neurogen/CAGE_PDBrainMap/output_dir/UW616/accepted_hits.plus.bw
CAGErev=~/neurogen/CAGE_PDBrainMap/output_dir/UW616/accepted_hits.minus.bw
RNAseq=~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.accepted_hits.normalized2.bw

# get mean RPM for promoter regions
awk '{OFS="\t"; if($4~/protein_coding/) {split($4,a,"___"); print a[3];}}' $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 > tmp
fgrep -wf tmp ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.isoforms.fpkm_tracking | cut -f1 > promoter.list

fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.promoter.tab
fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.promoter.tab
# RNAseq rpm at promoter regions
fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SNDA.promoter.tab
# RNAseq FPKM for Tx
fgrep -wf promoter.list ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SNDA_2_rep1.uniq.isoforms.fpkm_tracking | cut -f1,10 > RNAseq.SNDA.Tx.tab

# get mean RPM for permissive enhancers
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.cageEnhancer.tab
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.cageEnhancer.tab
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $RNAseq  stdin /dev/null -minMax -bedOut=RNAseq.SNDA.cageEnhancer.tab

# get mean RPM for HiTNE regions
bigWigAverageOverBed $CAGEfwd ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=CAGE.SNDA.fwd.HiTNE.tab
bigWigAverageOverBed $CAGErev ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=CAGE.SNDA.rev.HiTNE.tab
bigWigAverageOverBed $RNAseq  ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=RNAseq.SNDA.HiTNE.tab

# get mean RPM for eHiTNE regions
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SNDA.fwd.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SNDA.rev.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SNDA.eHiTNE.tab


CAGEfwd=~/eRNAseq/../CAGE/CAGE.fwd.bigwig
CAGErev=~/eRNAseq/../CAGE/CAGE.rev.bigwig
RNAseq=~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.accepted_hits.normalized2.bw

awk '{OFS="\t"; if($4~/protein_coding/) {split($4,a,"___"); print a[3];}}' $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 > tmp
fgrep -wf tmp ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.isoforms.fpkm_tracking | cut -f1 > promoter.list

fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.promoter.tab
fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.promoter.tab
# RNAseq rpm at promoter regions
fgrep -f promoter.list $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; if($4~/\.protein_coding/) {tss=($6=="+")?$2:$3; split($4,a,"___"); print $1,tss,tss+1, a[3];}}' | bedtools slop -b 500 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SN.promoter.tab
# RNAseq FPKM for Tx
fgrep -wf promoter.list ~/neurogen/rnaseq_PD/for_display/HC_UWA616_SN_6_rep1.amplified.uniq.isoforms.fpkm_tracking | cut -f1,10 > RNAseq.SN.Tx.tab

# get mean RPM for permissive enhancers
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.cageEnhancer.tab
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.cageEnhancer.tab
grep -v track ~/eRNAseq/../CAGE/permissive_enhancers.bed | cut -f1-4 | bigWigAverageOverBed $RNAseq  stdin /dev/null -minMax -bedOut=RNAseq.SN.cageEnhancer.tab

# get mean RPM for HiTNE regions
bigWigAverageOverBed $CAGEfwd ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=CAGE.SN.fwd.HiTNE.tab
bigWigAverageOverBed $CAGErev ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=CAGE.SN.rev.HiTNE.tab
bigWigAverageOverBed $RNAseq  ~/eRNAseq/eRNA.bed /dev/null -minMax -bedOut=RNAseq.SN.HiTNE.tab

# get mean RPM for eHiTNE regions
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $CAGEfwd stdin /dev/null -minMax -bedOut=CAGE.SN.fwd.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $CAGErev stdin /dev/null -minMax -bedOut=CAGE.SN.rev.eHiTNE.tab
awk '{OFS="\t"; if($1~/chr/ && ($6>=5 || $7>0 || $8>0 || $9>0 || $13>0 || $16>0)) {split($1,a,"_"); print a[1],a[2],a[3],$1}}' ~/eRNAseq/eRNA.characterize.xls | bigWigAverageOverBed $RNAseq stdin /dev/null -minMax -bedOut=RNAseq.SN.eHiTNE.tab

## merge together


join -1 1 -2 1 <(cat CAGE.SNDA.fwd.promoter.tab CAGE.SNDA.rev.promoter.tab | cut -f4,7 | sort) <(sort RNAseq.SNDA.Tx.tab) | awk '{OFS="\t"; print $1,($2>0)?$2:(-$2),$3,"SNDA","promoter"}'> CAGE.RNAseq.tab
#join -1 1 -2 1 <(cat CAGE.SNDA.fwd.promoter.tab CAGE.SNDA.rev.promoter.tab | cut -f4,7 | sort) <(cut -f4,7 RNAseq.SNDA.promoter.tab | sort) | awk '{OFS="\t"; print $1,($2>0)?$2:(-$2),$3,"SNDA","promoter"}'> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.cageEnhancer.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.cageEnhancer.tab | sort) | awk '{print $1,(($2-$3)>0)?$2:(-$3);}') <(cut -f4,7 RNAseq.SNDA.cageEnhancer.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","cageEnhancer"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.HiTNE.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.HiTNE.tab | sort) | awk '{print $1,(($2-$3)>0)?$2:(-$3);}') <(cut -f4,7 RNAseq.SNDA.HiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","HiTNE"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SNDA.fwd.eHiTNE.tab | sort) <(cut -f4,7 CAGE.SNDA.rev.eHiTNE.tab | sort) | awk '{print $1,(($2-$3)>0)?$2:(-$3);}') <(cut -f4,7 RNAseq.SNDA.eHiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SNDA","eHiTNE"}' >> CAGE.RNAseq.tab

#join -1 1 -2 1 <(cat CAGE.SN.fwd.promoter.tab CAGE.SN.rev.promoter.tab | cut -f4,7 | sort) <(cut -f4,7 RNAseq.SN.promoter.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","promoter"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(cat CAGE.SN.fwd.promoter.tab CAGE.SN.rev.promoter.tab | cut -f4,7 | sort) <(sort RNAseq.SN.Tx.tab) | awk '{OFS="\t"; print $1,$2,$3,"SN","promoter"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.cageEnhancer.tab | sort) <(cut -f4,7 CAGE.SN.rev.cageEnhancer.tab | sort) | awk '{print $1,($2>$3)?$2:$3}') <(cut -f4,7 RNAseq.SN.cageEnhancer.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","cageEnhancer"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.HiTNE.tab | sort) <(cut -f4,7 CAGE.SN.rev.HiTNE.tab | sort) | awk '{print $1,($2>$3)?$2:$3}') <(cut -f4,7 RNAseq.SN.HiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","HiTNE"}' >> CAGE.RNAseq.tab
join -1 1 -2 1 <(join -1 1 -2 1 <(cut -f4,7 CAGE.SN.fwd.eHiTNE.tab | sort) <(cut -f4,7 CAGE.SN.rev.eHiTNE.tab | sort) | awk '{print $1,($2>$3)?$2:$3}') <(cut -f4,7 RNAseq.SN.eHiTNE.tab | sort) | awk '{OFS="\t"; print $1,$2,$3,"SN","eHiTNE"}' >> CAGE.RNAseq.tab

##
R
df=read.table("CAGE.RNAseq.tab", header=F)
colnames(df)=c("id","CAGE","RNAseq","tissue_regio_for_CAGE","annotation_regions")
#dim(subset(df, CAGE==0 & RNAseq==0))
df=subset(df, CAGE>0 & RNAseq>0)
par(mfrow=c(2,4))
for(i in c("SN","SNDA")){
  for(j in c("promoter","cageEnhancer","HiTNE","eHiTNE"))
    plot(CAGE~RNAseq, data=subset(df, tissue_regio_for_CAGE==i & annotation_regions==j), log='xy', main=paste(i,j), pch=19, cex=0.6, col=ifelse(i=="SN",'#0000ff66','#228b2266'))
}

par(mfrow=c(1,4))
for(j in c("promoter","cageEnhancer","HiTNE","eHiTNE")){
  plot(subset(df, tissue_regio_for_CAGE=='SN' & annotation_regions==j, select = CAGE)[,1], subset(df, tissue_regio_for_CAGE=='SNDA' & annotation_regions==j, select = CAGE)[,1], xlab="FANTOM CAGE (SN)", ylab="BRAINCODE CAGE (SNDA)", log='xy', main=j, pch=19, cex=0.6, col='#00000066')
}
    



plot(CAGE~RNAseq, data=subset(df, tissue_regio_for_CAGE=="SN" & annotation_regions=="promoter"), log='xy', pch=19, cex=0.6, col='#0000ff66')

require(ggplot2)
p <- ggplot(df, aes(CAGE, RNAseq)) + geom_point() + scale_x_log10() + scale_x_log10() 
p + facet_grid(tissue_regio_for_CAGE ~ annotation_regions, scales='free')
