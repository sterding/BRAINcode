# ===========================================================================
# output:
# 1. blacklist.bed: The background region to measure transcriptional noise, defined as genomic regions excluding the known regions with RNA activities (known exons+/-500bp, rRNA, CAGE-defined enhancers, promoters)
# 2. toExclude.bed: the above blacklist - known enhancers regions (e.g. FANTOM5 permissive enhancers)
# ===========================================================================

cd ~/eRNAseq

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

## deadregion: random regions without RNAseq and CAGE
# ----------------------------------------------------

complementBed -i ~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bedGraph -g <(sort -k1,1 $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size) > noRNAseq.HCILB_SNDA.bed &
complementBed -i ~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HC_PY.bedGraph -g <(sort -k1,1 $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size) > noRNAseq.HC_PY.bed &
complementBed -i ~/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.HC_nonNeuron.bedGraph -g <(sort -k1,1 $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size) > noRNAseq.HC_nonNeuron.bed &

bigWigToBedGraph ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.fwd.bigwig stdout | bedtools slop -r 100 -l 0 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | mergeBed > ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.fwd.bigwig.r100bp.bed
bigWigToBedGraph ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.rev.bigwig stdout | bedtools slop -l 100 -r 0 -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | mergeBed > ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.rev.bigwig.l100bp.bed

intersectBed -a noRNAseq.HCILB_SNDA.bed -b noRNAseq.HC_PY.bed | intersectBed -a - -b noRNAseq.HC_nonNeuron.bed | intersectBed -a - -b <(cat ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.fwd.bigwig.r100bp.bed ~/eRNAseq/externalData/CAGE/CAGE.FANTOM5.total.rev.bigwig.l100bp.bed | sortBed | mergeBed) -v | intersectBed -a - -b $GENOME/Annotation/Variation/rmsk.hg19.bed -v | intersectBed -a - -b $GENOME/Annotation/Genes/exons.meta.bed -v | awk '($3-$2)>150' > noRNAseq.noCAGE.noRepeat.noExons.bed

cd -