# ===========================================================================
#: background region to measure transcriptional noise: genomic regions excluding the known regions with RNA activities (known exons+/-500bp, rRNA, CAGE-defined enhancers, promoters)
# ===========================================================================

cd ~/projects/PD/results/eRNA

ANNOTATION=$GENOME/Annotation/Genes

cat $ANNOTATION/gencode.v19.annotation.bed12 $ANNOTATION/knownGene.bed12 $ANNOTATION/NONCODEv4u1_human_lncRNA.bed12 | bed12ToBed6 | cut -f1-3 | grep -v "_" | slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 > /tmp/bg.bed
cut -f1-3 $ANNOTATION/rRNA.bed >> /tmp/bg.bed  # rRNA
# +/-500bp flanking around the CAGE-predicted TSS (downloaded from: http://fantom.gsc.riken.jp/5/datafiles/latest/extra/TSS_classifier/)
grep -v track ~/projects/PD/results/eRNA/externalData/CAGE/TSS_human.bed | grep -v "211,211,211" | cut -f1-3 | grep -v "_" | slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 >> /tmp/bg.bed 
#cat $ANNOTATION/SINE.bed $ANNOTATION/LINE.bed | cut -f1-3 >> /tmp/bg.bed  # LINE and SINE
cat $ANNOTATION/hg19.gap.bed >> /tmp/bg.bed  # genomic gap
cat /tmp/bg.bed | sortBed | mergeBed -i - > toExclude.bed
grep -v track ~/projects/PD/results/eRNA/externalData/CAGE/permissive_enhancers.bed | cut -f1-3 >> /tmp/bg.bed # CAGE-enhancer
cat /tmp/bg.bed | sortBed | mergeBed -i - > blacklist.bed

cd -