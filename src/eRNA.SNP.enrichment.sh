# Usage: 
# for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 
# for i in HCILB_SNDA HC_TCPY HC_MCPY HC_FB HC_PBMC; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 

# alternatively,
## for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh PLINK $i; done 
# bsub -n 1 -q normal -J HCILB_SNDA bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP HCILB_SNDA

type=$1
samplegroup=$2

# debug
# type='SNAP'; samplegroup='HCILB_SNDA'

cd ~/eRNAseq/$samplegroup

## pre-steps to get LD data for GWAS SNPs: 
## Ref ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/README.txt

[ "$type" == "SNAP" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed  # in the final figure we used SNAP
#[ "$type" == "PLINK" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.PLINK.LD_w250.r2_0.8.bed

## extract all autosomal.associations
[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sort -u > $snps_in_LD.autosomal.associations.bed

# number of gwas SNPs
wc -l $snps_in_LD
# number of diseases/traits
cut -f7 $snps_in_LD |sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort -u | wc -l
# number of associations
wc -l $snps_in_LD.autosomal.associations.bed

## random regions
#bedtools shuffle -excl <(cat blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -noOverlapping -i eRNA.bed -g $GENOME/Annotation/Genes/ChromInfo.txt > eRNA.random.bed
bedtools random -n 100000 -l 400 -seed 1234 -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(cat ../blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -v > eRNA.random.bed
# wc -l eRNA.random.bed --> N=72202 (close to 71469 lines of eRNA.bed)

## mRNA inner exons
grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed | awk '{if(id!=$4) id=$4; else print}' | sort -k4,4 -k2,2nr | awk '{if(id!=$4) id=$4; else print}' > $GENOME/Annotation/Genes/mRNA.innner.exon.bed

## promoter -- [-200,+200] of protein-coding GENCODE v19 TSS
grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400}' > $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed

echo "## overlapped SNPs with each dataset"
### ##################
echo "# all"
cat $snps_in_LD.autosomal.associations.bed | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.all
echo "# eRNA"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE
echo "# eRNA - class I"
[ -e eRNA.classI.bed ] || awk 'NR>1{OFS="\t"; split($1,a,"_"); if($28==1) print a[1],a[2],a[3],$1}' eRNA.characterize.xls > eRNA.classI.bed
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.classI.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE1
echo "# eRNA - class II"
[ -e eRNA.classII.bed ] || awk 'NR>1{OFS="\t"; split($1,a,"_"); if($28==2) print a[1],a[2],a[3],$1}' eRNA.characterize.xls > eRNA.classII.bed
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.classII.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE2
echo "# eRNA - class III"
[ -e eRNA.classIII.bed ] || awk 'NR>1{OFS="\t"; split($1,a,"_"); if($28==3) print a[1],a[2],a[3],$1}' eRNA.characterize.xls > eRNA.classIII.bed
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.classIII.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE3
echo "# eRNA - class I+II"
[ -e eRNA.classInII.bed ] || awk 'NR>1{OFS="\t"; split($1,a,"_"); if($28<3) print a[1],a[2],a[3],$1}' eRNA.characterize.xls > eRNA.classInII.bed
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.classInII.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE1n2
echo "# eRNA-private"
[ -e eRNA.private.major.bed ] && intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.private.major.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE-private
#[ -e eRNA.private.minor.bed ] && intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.private.minor.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.private.minor.HTNE
echo "# mRNA inner exons"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $GENOME/Annotation/Genes/mRNA.innner.exon.bed  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.exon
echo "# promoter -- [-200,+200] of protein-coding GENCODE v19 TSS"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.promoter
echo "# randomly sampling"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.random.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.random

## add enhancers defined by other features (12/08/2017)
EXTERNAL_FEATURE=~/eRNAseq/externalData
echo "# chromHMM-brain"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/Segment/15_coreMarks_segments.E7enhancer.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.chromHMM_brain
echo "# chromHMM-cellline"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/Segment/wgEncodeBroadHmm.strongEnhancer.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.chromHMM_cellline
echo "# DNase"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/DNase/regions_enh_merged.brain.narrowPeak -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.DNase
echo "# CAGE"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/CAGE/permissive_enhancers.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.CAGE
echo "# TFBS hotspot"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.hotspot.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.TFhotspot
echo "# P300"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b <(awk '$4=="EP300"' $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed) -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.P300
echo "# Conservation"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $EXTERNAL_FEATURE/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HCNE


#rm $snps_in_LD.autosomal.associations.bed

echo "## total overlapped SNPs count with each dataset"
### ##################
echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.$type.counts.summary
echo "HTNE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HTNE1n2" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.classInII.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HTNE1" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.classI.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HTNE2" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.classII.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HTNE3" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.classIII.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HTNE-private" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.private.major.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "exon" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $GENOME/Annotation/Genes/mRNA.innner.exon.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "promoter" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "random" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.random.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
## add enhancers defined by other features (12/08/2017)
echo "chromHMM_cellline" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/Segment/wgEncodeBroadHmm.strongEnhancer.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "chromHMM_brain" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/Segment/15_coreMarks_segments.E7enhancer.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "DNase" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/DNase/regions_enh_merged.brain.narrowPeak -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "CAGE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/CAGE/permissive_enhancers.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "TFhotspot" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.hotspot.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "P300" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b <(awk '$4=="EP300"' $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed) -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "HCNE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $EXTERNAL_FEATURE/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary

#echo "## Fisher test and make plot"  # move out of the script now
### ##################
#Rscript $pipeline_path/src/eRNA.SNP.enrichment.R $type $samplegroup