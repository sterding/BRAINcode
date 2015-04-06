## script to generate data to draw aggregation plot for enhancer regions 

ANNOTATION=$GENOME/Annotation/Genes

cd ~/projects/PD/results/eRNA/externalData/

# regions to exclude
#=====================
# regions to exclude: 1) 500nt flanking of GENCODE exons, 2) 500nt flanking of UCSC KnowGene exons (since some UCSC genes are not annotated in GENCODE, e.g. chr1:665045-665121), 3) rRNA from repeatMasker
cat $ANNOTATION/gencode.v19.annotation.bed12 $ANNOTATION/knownGene.bed12 | bed12ToBed6 | cut -f1-3| grep -v "_" |slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 500 | cat - <(cut -f1-3 $ANNOTATION/rRNA.bed) | sortBed | mergeBed > toExclude.bed

#============================================================
# Enhancers defined by FANTOM5 CAGE
#============================================================
cd ~/projects/PD/results/eRNA/externalData/CAGE
curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed | awk '{OFS="\t"; peak=$8; if(peak>500) print $1, peak-500, peak+500;}' > CAGE.distal.bed

wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.cell_line.hCAGE/epitheloid%20carcinoma%20cell%20line:%20HelaS3%20ENCODE,%20biol_rep1.CNhs12325.10815-111B5.hg19.ctss.fwd.bw
wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.cell_line.hCAGE/epitheloid%20carcinoma%20cell%20line:%20HelaS3%20ENCODE,%20biol_rep1.CNhs12325.10815-111B5.hg19.ctss.rev.bw
wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.tissue.hCAGE/substantia%20nigra,%20adult,%20donor10252.CNhs12318.10158-103A5.hg19.ctss.fwd.bw
wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.tissue.hCAGE/substantia%20nigra,%20adult,%20donor10252.CNhs12318.10158-103A5.hg19.ctss.rev.bw

ln -fs substantia*ctss.fwd.bw ../CAGE.fwd.bigwig
ln -fs substantia*ctss.rev.bw ../CAGE.rev.bigwig

#============================================================
# Enhancers defined by ENCODE histone modifications
# e.g. active enhancers = H3K4me1 + H3k27ac - H3K27me3
# silent enhancers = H3K4me1 - H3k27ac + H3K27me3
#============================================================
cd ~/projects/PD/results/eRNA/externalData/ENCODE
url=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone
for i in K562 Gm12878 H1hesc;
do
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz;
    ##wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k04me1StdPk.broadPeak.gz;
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz;
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz;
     
    intersectBed -a <(zcat wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz) -b <(zcat wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz) -wa -u | intersectBed -a - -b <(zcat wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz) -v -wa | awk -vi=$i '{OFS="\t"; $4=i; print}' > $i.silent.enhancer.broadPeak
    intersectBed -a <(zcat wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz) -b <(zcat wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz) -wa -u | intersectBed -a - -b <(zcat wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz) -v -wa | awk -vi=$i '{OFS="\t"; $4=i; print}' > $i.active.enhancer.broadPeak
    echo $i;
done

# P300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsSydhK562Gata1UcdUniPk.narrowPeak.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsBroadK562P300UniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "K562";}' | intersectBed -a - -b K562.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.bed -v -wa > K562.active.enhancer.broadPeak.distal.wP300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibGm12878P300Pcr1xUniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "Gm12878";}' | intersectBed -a - -b Gm12878.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.bed -v -wa > Gm12878.active.enhancer.broadPeak.distal.wP300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibH1hescP300V0416102UniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "H1hesc";}' | intersectBed -a - -b H1hesc.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.bed -v -wa > H1hesc.active.enhancer.broadPeak.distal.wP300

# merge
sort -k1,1 -k2,2n *active.enhancer.broadPeak.distal.wP300 > ENCODE.distal.bed

# bigwig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k04me1StdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k27acStdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k27me3StdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k4me3StdSig.bigWig

ln -s wgEncodeBroadHistoneHelas3H3k04me1StdSig.bigWig ENCODE.H3k4me1.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k27acStdSig.bigWig ENCODE.H3k27ac.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k27me3StdSig.bigWig ENCODE.H3k27me3.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k4me3StdSig.bigWig ENCODE.H3k4me3.bigwig

# use Brain_Substantia_Nigra data
# download from: http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?search=substantia+nigra&display=50

wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM772nnn/GSM772898/suppl/GSM772898_BI.Brain_Substantia_Nigra.H3K4me1.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM772898_BI.Brain_Substantia_Nigra.H3K4me1.149.bw
ln -s !$ Histone.SN.H3K4me1.bigwig
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM772nnn/GSM772901/suppl/GSM772901_BI.Brain_Substantia_Nigra.H3K4me3.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM772901_BI.Brain_Substantia_Nigra.H3K4me3.149.bw
ln -s !$ Histone.SN.H3K4me3.bigwig
curl -s ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM772nnn/GSM772937/suppl/GSM772937_BI.Brain_Substantia_Nigra.H3K27me3.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM772937_BI.Brain_Substantia_Nigra.H3K27me3.149.bw
ln -s !$ Histone.SN.H3K27me3.bigwig
curl -s ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1112nnn/GSM1112778/suppl/GSM1112778_BI.Brain_Substantia_Nigra.H3K27ac.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM1112778_BI.Brain_Substantia_Nigra.H3K27ac.149.bw
ln -s !$ Histone.SN.H3K27ac.bigwig
curl -s ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM916nnn/GSM916015/suppl/GSM916015_BI.Brain_Substantia_Nigra.H3K36me3.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM916015_BI.Brain_Substantia_Nigra.H3K36me3.149.bw
ln -s !$ Histone.SN.H3K36me3.bigwig
curl -s ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM669nnn/GSM669977/suppl/GSM669977_BI.Brain_Substantia_Nigra.H3K9ac.112.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM669977_BI.Brain_Substantia_Nigra.H3K9ac.112.bw
ln -s !$ Histone.SN.H3K9ac.bigwig

#============================================================
# Enhancers defined by TFBS HOT region
#============================================================
cd ../TFBS

# clustered TFBS (count all Peaks for 161 transcription factors in 91 cell types)
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip -c wgEncodeRegTfbsClusteredV3.bed.gz | awk '{OFS="\t"; print $1,$2,$3,$4,0,"+",$2,$3,"0",$6,$7,$8}' > wgEncodeRegTfbsClusteredV3.bed12
# per TF per cell type
#sort -k1,1 wgEncodeRegTfbsClusteredV3.bed12 | bedItemOverlapCount hg19 -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > wgEncodeRegTfbsClusteredV3.bed12.bg
# only per TF (e.g. if a TF occur in >1 cell types, it's only counted once)
sort -k1,1 wgEncodeRegTfbsClusteredV3.bed12 | bedItemOverlapCount hg19 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > wgEncodeRegTfbsClusteredV3.bg
bedGraphToBigWig wgEncodeRegTfbsClusteredV3.bg $ANNOTATION/ChromInfo.txt wgEncodeRegTfbsClusteredV3.bw
ln -fs wgEncodeRegTfbsClusteredV3.bw TFBS.bigwig

### Option1: use the hot region from site below:
#curl -s http://stanford.edu/~claraya/metrn/data/hot/regions/hs/maphot_hs_selection_reg_cx_occP05_any.bed | cut -f1-6  > maphot_hs_selection_reg_cx_occP05_any.bed6
#intersectBed -a maphot_hs_selection_reg_cx_occP05_any.bed6 -b ../toExclude.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > TFBS.distal.bed

## Option2: use the ENCODE TFBS cluster data, any region with >5 TFs (including P300) in any cell lines
awk '{OFS="\t"; if($4>5) print $1,$2,$3,".",$4}' wgEncodeRegTfbsClusteredV3.bg | mergeBed -scores max | intersectBed -a stdin -b <( grep P300 wgEncodeRegTfbsClusteredV3.bed12 | cut -f1-3) -u > TFBS.distal.bed

# length check up
awk '{print $3-$2}' TFBS.distal.bed | textHistogram stdin -binSize=50 -maxBinCount=50

#============================================================
# Enhancers defined by conservation
#============================================================
cd ../Conservation
# Ancora HCNEs
curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_hg19_danRer7_70pc_50col.bed.gz | gunzip | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp
ln -sf HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp Conservation.distal.bed


## bigwig file
#curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_density_hg19_danRer7_70pc_50col.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size HCNE_density_hg19_danRer7_70pc_50col.bw
#ln -s HCNE_density_hg19_danRer7_70pc_50col.bw HCNE.distal.bigwig

# phyloP
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate ./
zcat vertebrate/chr*.gz | gzip -c > vertebrate/phyloP46way.wigFix.gz
wigToBigWig vertebrate/phyloP46way.wigFix.gz $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size vertebrate/phyloP46way.wigFix.bigwig

ln -sf vertebrate/phyloP46way.wigFix.bigwig Conservation.bigwig

## 29 mammals conservation (in hg18)
#curl -s http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz | gunzip > 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt
#intersectBed -a way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt -b ../toExclude.bed -v -wa
#
## bigwig file
#wget http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/omega.12mers.wig.gz
#wigToBigWig omega.12mers.wig $GENOME/Sequence/WholeGenomeFasta/hg18.chrom.size omega.12mers.hg18.bigwig

#============================================================
# Enhancers defined by DNase
#============================================================
cd ../DNase
# ref: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=377332155_oXR9pqyaZLzzFHxgO4t0YzyQseGN&c=chr1&g=wgEncodeRegDnaseClusteredV2
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV2.bed.gz
zcat wgEncodeRegDnaseClusteredV2.bed.gz | awk '$4>20 && $5>800' | intersectBed -a - -b ../toExclude.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > DNase.distal.bed  # at least occuring in 20 cell lines (out of 125), and 800 cluster score (out of 1000) --> 100180 regions remained

#wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseH1hescPkRep1.narrowPeak.gz
#zcat wgEncodeUwDnaseH1hescPkRep1.narrowPeak.gz | intersectBed -a - -b ../toExclude.bed -v | cut -f1-3 > DNase.Helas3.bed

wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseHelas3RawRep1.bigWig
ln -s wgEncodeUwDnaseHelas3RawRep1.bigWig DNase.bigwig

curl -s http://www.genboree.org/EdaccData/Current-Release/sample-experiment/Fetal_Brain/Chromatin_Accessibility/UW.Fetal_Brain.ChromatinAccessibility.H-24510.DNase.DS20780.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size UW.Fetal_Brain.ChromatinAccessibility.H-24510.DNase.DS20780.bw
ln -fs UW.Fetal_Brain.ChromatinAccessibility.H-24510.DNase.DS20780.bw DNase.bigwig

#============================================================
# Enhancers defined by RNAseq
#============================================================
# uniq mapping RNAseq with (1) RPM>0.5 (2) non-genic and (3) length>=50   --> uniq/accepted_hits.normalized.eRNA.bed in each sample folder
cd ../RNAseq

awk '$4>0.05' /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bedGraph | mergeBed -d 100 | intersectBed -a - -b ../toExclude.bed -v | awk '($3-$2)>=200' > HC_SNDA.trimmedmean.uniq.normalized.eRNA.strong.bed
awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' HC_SNDA.trimmedmean.uniq.normalized.eRNA.strong.bed > RNAseq.distal.bed

ln -s HC_SNDA.trimmedmean.uniq.normalized.bw RNAseq.distal.bigwig

#============================================================
# defined by Methylation
#============================================================
cd ../Methylation
curl -s ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM916nnn/GSM916086/suppl/GSM916086_BI.Brain_Substantia_Nigra.RRBS.149.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size GSM916086_BI.Brain_Substantia_Nigra.RRBS.149.bw

ln -fs GSM916086_BI.Brain_Substantia_Nigra.RRBS.149.bw Methylation.bigwig


#============================================================
# get bin signal
#============================================================
cd ~/neurogen/rnaseq_PD/results/eRNA/externalData

for i in Histone CAGE TFBS Conservation DNase RNAseq;
do
    for j in */*bigwig; do
        echo $i, $j;
        [ -f $j.$i.bins ] || toBinRegionsOnBigwig.sh $j $i/$i.distal.bed 200 > $j.$i.bins &
    done
    #for expressed regions only (for DoD grant)
    [ -f $i/${i}_n_RNAseq.distal.bed ] || intersectBed -a $i/$i.distal.bed -b RNAseq/RNAseq.distal.bed -u > $i/${i}_n_RNAseq.distal.bed
    [ -f RNAseq/RNAseq.${i}_n_RNAseq.bins ] || toBinRegionsOnBigwig.sh RNAseq/RNAseq.bigwig $i/${i}_n_RNAseq.distal.bed 200 > RNAseq/RNAseq.${i}_n_RNAseq.bins &
done


#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/neurogen/pipeline/RNAseq/src/enhancer_analysis.aggPlot.R

#============================================================
# draw venn diagram
#============================================================
RNAseq="RNAseq/RNAseq.distal.bed"
CAGE="CAGE/CAGE.distal.bed"
Histone="Histone/Histone.distal.bed"

wc -l $RNAseq  #46157
wc -l $CAGE  # 43011
wc -l $Histone  #6009

intersectBed -a $RNAseq -b $CAGE -f 0.5 -r | wc -l  # 863
intersectBed -a $RNAseq -b $Histone -f 0.5 -r | wc -l  # 132
intersectBed -a $CAGE -b $Histone -f 0.5 -r | wc -l  # 1496

intersectBed -a $RNAseq -b $CAGE -wa | wc -l  # 1569
intersectBed -a $RNAseq -b $Histone -wa | wc -l  # 224
intersectBed -a $CAGE -b $Histone -wa | wc -l  # 1933
intersectBed -a $CAGE -b $Histone -wa | intersectBed -a - -b $RNAseq -wa | wc -l  # 49

R
require("venneuler")
pdf("vendiagram.CAGE.RNAseq.Hist.pdf")
plot(venneuler(c(RNAseq=46157/46157, CAGE=43011/46157, Histone=6009/46157, "RNAseq&CAGE"=1569/46157, "RNAseq&Histone"=224/46157, "CAGE&Histone"=1933/46157 ,"RNAseq&CAGE&Histone"=490/46157)))
dev.off()