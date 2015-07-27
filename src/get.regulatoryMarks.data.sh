#============================================================
# Script to download bigwig files for various marks
# Usage: bash $pipeline_path/src/get.regulatoryMarks.data.sh
#============================================================

ANNOTATION=$GENOME/Annotation/Genes
externalData=~/projects/PD/results/eRNA/externalData
[ -d $externalData ] || mkdir $externalData
cd $externalData

#========================================
# get bigwig and enhancer for characteristic marks
#========================================

echo "getting CAGE..."
# ----------------------------------
[ -d $externalData/CAGE ] || mkdir $externalData/CAGE
cd $externalData/CAGE

curl -s http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.tissue.hCAGE/substantia%20nigra,%20adult,%20donor10252.CNhs12318.10158-103A5.hg19.ctss.fwd.bw > CAGE.FANTOM5.SN.fwd.bigwig
curl -s http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.tissue.hCAGE/substantia%20nigra,%20adult,%20donor10252.CNhs12318.10158-103A5.hg19.ctss.rev.bw > CAGE.FANTOM5.SN.rev.bigwig
curl -s http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.fwd.bw > CAGE.FANTOM5.total.fwd.bigwig
curl -s http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.rev.bw > CAGE.FANTOM5.total.rev.bigwig

# permissive enhancers
curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed > permissive_enhancers.bed
# permissive TSS
curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/TSS_classifier/TSS_human.bed.gz | gunzip > TSS_human.bed

echo "getting Histone..."
# ---------------------------------
[ -d $externalData/Histone ] || mkdir $externalData/Histone
cd $externalData/Histone

curl -s http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K4me1.pval.signal.bigwig > H3K4me1.Roadmap.SN.pval.bigwig
curl -s http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K27ac.pval.signal.bigwig > H3K27ac.Roadmap.SN.pval.bigwig
curl -s http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K9ac.pval.signal.bigwig  > H3K9ac.Roadmap.SN.pval.bigwig


echo "getting TFBS HOT region"
# ---------------------------------
[ -d $externalData/TFBS ] || mkdir $externalData/TFBS
cd $externalData/TFBS

# clustered TFBS (count all Peaks for 161 transcription factors in 91 cell types)
# only per TF (e.g. if a TF occur in >1 cell types, it's only counted once)

curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz | gunzip | awk '{OFS="\t"; print $1,$2,$3,$4,0,"+",$2,$3,"0",$6,$7,$8}' > wgEncodeRegTfbsClusteredV3.bed12
sort -k1,1 wgEncodeRegTfbsClusteredV3.bed12 | bedItemOverlapCount hg19 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > wgEncodeRegTfbsClusteredV3.bg
bedGraphToBigWig wgEncodeRegTfbsClusteredV3.bg $ANNOTATION/ChromInfo.txt TFBS.ENCODE.all.count.bigwig


echo "getting Conservation"
# ---------------------------------
[ -d $externalData/Conservation ] || mkdir $externalData/Conservation
cd $externalData/Conservation

# Ancora HCNEs
curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_hg19_danRer7_70pc_50col.bed.gz | gunzip | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp
ln -sf HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp Conservation.HCNE.distal.bed

# # phyloP
# # http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/
# rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate ./
# zcat vertebrate/chr*.gz | gzip -c > vertebrate/phyloP46way.wigFix.gz
# rm vertebrate/chr*.gz
# wigToBigWig vertebrate/phyloP46way.wigFix.gz $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size Conservation.phyloP46way.bigwig


echo "getting DNase"
# ---------------------------------
[ -d $externalData/DNase ] || mkdir $externalData/DNase
cd $externalData/DNase

# # ref: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=377332155_oXR9pqyaZLzzFHxgO4t0YzyQseGN&c=chr1&g=wgEncodeRegDnaseClusteredV2
# # at least occuring in 20 cell lines (out of 125), and 800 cluster score (out of 1000) --> 100180 regions remained
# curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV2.bed.gz | gunzip | awk '$4>20 && $5>800' | intersectBed -a - -b ../toExclude.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > DNase.distal.bed 

# ENCODE Hela
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseHelas3RawRep1.bigWig > DNase.ENCODE.Hela.bigwig

# fetal brain male
curl -s http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-DNase.pval.signal.bigwig > DNase.Roadmap.fBrain.pval.bigwig
curl -s http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E081-DNase.macs2.narrowPeak.gz | gunzip > DNase.Roadmap.fBrain.narrowPeak

echo "getting RNAseq"
# ---------------------------------
[ -d $externalData/RNAseq ] || mkdir $externalData/RNAseq
cd $externalData/RNAseq
ln -s /data/neurogen/rnaseq_PD/results2/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bw RNAseq.BRAINCODE.HCILB_SNDA.bigwig