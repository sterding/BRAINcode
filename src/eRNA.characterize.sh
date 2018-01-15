## script to characterize eRNAs with the following features
# $pipeline_path/src/eRNA.characterize.sh HCILB_SNDA
#
#Features	            Data type	    Description
#========               =========       ==================================
#dis2TSS	            integer	        distance between the middle points of HiTNE and the nearest TSS, in bp. If HiTNE is located intronic, it's plus; otherwise it's minus. 
#RPKM	                float	        normalized expression level, calculated in the same way as RPKM.
#RPM	                float	        reads density at the summit position, normalized to total mapped reads in million. 
#readsCount	            integer	        raw reads count mapped to the HiTNE
#normCpG	            float	        normalized CpG score
#nTFBS	                integer	        number of distinct TFs bound to the HiTNE, based on ENCODE ChIPseq cluster (wgEncodeRegTfbsClusteredV3)
#P300	                boolean	        if the P300 binding site found in the HiTNE
#enhancer_CAGE	        boolean	        if overlap with any CAGE-defined permissive enhancers
#enhancer_histone	    boolean	        if overlap with any histone marks-defined enhancers (chromHMM states of E6|E7|E12 from substantial nigro)
#enhancer_VISTA	        boolean	        if overlap with any tested enhancers (positive enhancers from VISTA enhancer database)
#DNase	                float	        maximal DNase density of fetal brain from Roadmap Epigenomics
#bDNase	                boolean	        if overlap with DNase cluster from ENCODE (wgEncodeRegDnaseClustered V2)
#DNaseROADMAP           boolean         if overlap with DNase narrow peak from Roadmap brain samples 
#conservation	        float	        mean phastCons score for the HiTNE region
#bConserved2zf	        boolean	        if conserved to zebrafish (defined by existence of chain alignment between human and zebrafish)
#bHCNE	                boolean	        if overlapping with any HCNEs (HCNE_hg19_danRer7_70pc_50col from Ancora)
#GWAS	                integer	        number of GWAS SNPs in HiTNE
#bGWAS	                boolean	        if any GWAS SNPs in HiTNE
#eSNP	                integer	        number of eQTL SNPs in HiTNE
#beSNP	                boolean	        if any eQTL SNPs in HiTNE
#nHostgene	            integer	        number of HiTNEs in the host gene. 0 for intergenic HiTNE
#lenHostgene	        integer	        length of host gene. 0 for intergenic HiTNE
#enhancer_CAGE.HeLa	    boolean	        if overlap with any CAGE-defined enhancers in HeLa cell line  -- for in vitro validation purpose only
#enhancer_histone.HeLa	boolean	        if overlap with any histone marks-defined enhancers (chromHMM states of E6|E7|E12) from HeLa Cell -- for in vitro validation purpose only

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

SAMPLE_GROUP=$1
inputBG=/data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph
EXTERNAL_FEATURE=~/eRNAseq/externalData

cd ~/eRNAseq/$SAMPLE_GROUP

inputbed=eRNA.bed

[ -e eRNA.random.bed ] || bedtools random -n 100000 -l 400 -seed 1234 -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(cat ../blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -v > eRNA.random.bed

# ====================================
## dis2TSS (distance btw middle of HiTNE and the nearest TSS)
# ====================================
#fgrep -w gene $GENOME/Annotation/Genes/gencode.v19.annotation.gtf | sed 's/[;"]//g'  | awk '{OFS="\t"; print $1, $4-1, $5, $18"___"$10"___"$14, 0, $7}' > $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed
# have to deal with intronic and intergenic separately
# intronic ones (if located in two genes' intron, just randomly pick the first hit in the file.)
TMPFILE=`mktemp /tmp/example.XXXXXXXXXX`
awk '{OFS="\t"; mid=int(($3+$2)/2); print $1, mid, mid+1,$4}' $inputbed | closestBed -a - -b <(sortBed -i $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed) -d -t first | awk '$11==0' | awk '{OFS="\t"; tss=($10=="+")?$6:$7; d=tss-$2; if(d<0) d=-d; print $4, d;}' > $TMPFILE
# intergenic ones
awk '{OFS="\t"; mid=int(($3+$2)/2); print $1, mid, mid+1,$4}' $inputbed | closestBed -a - -b <(sortBed -i $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed) -d -t first | awk '$11!=0' | cut -f1-4 | sort -k1,1 -k2,2n -u | closestBed -a - -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed | sortBed) -D b -t first | awk '{OFS="\t"; print $4,($11<0)?$11:-$11;}' >> $TMPFILE

sort -k1,1 $TMPFILE > eRNA.f01.dis2TSS.txt

# ====================================
## RPKM (Note: this RPKM is different from the normal RPKM, they might have a factor of read length difference)
# ====================================
# mean0: average over bases with non-covered bases counting as zeroes
bigWigAverageOverBed ${inputBG/bedGraph/bw} $inputbed stdout | cut -f1,5 | sort -k1,1 | awk '{OFS="\t"; print $1, $2*1000+0}' > eRNA.f02.RPKM.txt

# ====================================
### RPM 
# ====================================
bigWigAverageOverBed ${inputBG/bedGraph/bw} $inputbed stdout -minMax | cut -f1,8 | sort -k1,1  > eRNA.f03.RPM.txt

# ====================================
# normalized CpG score
# ====================================
bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed $inputbed -fo eRNA.seq.tab
$pipeline_path/bin/getNormalizedCpGscore.awk eRNA.seq.tab | sort -k1,1  > eRNA.f05.CpG.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 eRNA.f05.CpG.tab

bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed eRNA.random.bed -fo random.seq.tab
$pipeline_path/bin/getNormalizedCpGscore.awk random.seq.tab | sort -k1,1 > random.f05.CpG.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 random.f05.CpG.tab

grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400,$4}' | bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo promoters.seq.txt
$pipeline_path/bin/getNormalizedCpGscore.awk promoters.seq.txt | sort -k1,1 > promoters.f05.CpG.txt
#textHistogram -col=2 -real -binSize=0.02 -maxBinCount=50 -minVal=0 promoters.f05.CpG.tab

# ====================================
# TFBS
# ====================================
## TFBS data for 161 transcription factors in 91 cell types was downloaded from ENCODE
#curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz | gunzip | awk '{OFS="\t"; $6=1+gsub(",",",",$6); print}' > ../TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed
# number of TF ChIP-seq peaks in the region (only if it occurs in at least one cell lines ## only TF peaks supported by >=5 cell lines are counted)
awk '$6>=1' $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $inputbed -b stdin -c | sort -k4,4 | cut -f4,5 > eRNA.f06.TFBS.txt

# Define TF hotspot (regions with >=5 different TFs bound)
sort -k1,1 -k2,2n $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | mergeBed -i stdin -c 4 -o count_distinct | awk '$4>=5' > $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.hotspot.bed 

# ====================================
# P300 binding
# ====================================
# if any P300 biding site found in the region
awk '$4=="EP300"' $EXTERNAL_FEATURE/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $inputbed -b stdin -c | sort -k4,4 | cut -f4,5 > eRNA.f07.P300.txt

# ====================================
# CAGE-defined enhancers
# ====================================
## download from FANTOM5 premissive enhancers:
# curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed > ../CAGE/permissive_enhancers.bed  #N=43011
# curl -s http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz | gunzip > ../CAGE/human_permissive_enhancers_phase_1_and_2.bed  # N=65423

intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/permissive_enhancers.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f08.CAGEenhancer.txt
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/human_permissive_enhancers_phase_1_and_2.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f08.CAGEenhancer.1n2.txt

# ====================================
# overlap with histone-defined enhancers
# ====================================

# Roadmap Epigenomics enhancers (http://egg2.wustl.edu/roadmap/web_portal/meta.html)
# download brain chromHMM segmentation data for 10 brain tissues:
#Brain Hippocampus Middle : 
#Brain Substantia Nigra
#Brain Anterior Caudate
#Brain Cingulate Gyrus
#Brain Inferior Temporal Lobe
#Brain Angular Gyrus
#Brain_Dorsolateral_Prefrontal_Cortex
#Brain Germinal Matrix
#Fetal Brain Female
#Fetal Brain Male
#for i in E067 E068 E069 E070 E071 E072 E073 E074 E082 E081; do echo $i; curl -s http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/${i}_15_coreMarks_segments.bed > ../Segment/${i}_15_coreMarks_segments.bed; done
# E6, E7, and E12 for enhancer (http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state)
# cat externalData/Segment/E*_15_coreMarks_segments.bed | awk '$4~/E6|E7|E12/' > externalData/Segment/15_coreMarks_segments.E6E7E12.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Segment/15_coreMarks_segments.E6E7E12.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f09.chromHMM_brain.txt

# strong enhancer (E7)
# cat externalData/Segment/E*_15_coreMarks_segments.bed | awk '$4=="E7"' > externalData/Segment/15_coreMarks_segments.E7enhancer.bed
# active promoter (E1)
# cat externalData/Segment/E*_15_coreMarks_segments.bed | awk '$4=="E1"' > externalData/Segment/15_coreMarks_segments.E1promoter.bed

# ====================================
# overlap with VISTA enhancers
# ====================================
intersectBed -a $inputbed -b <(cut -f1-4,7 $EXTERNAL_FEATURE/VISTA/hg19.tested_regions.bed) -wao  | sort -k4,4 | awk '{OFS="\t"; print $4, ($6==-1)?"NA":$8"|"$9}' | groupBy -g 1 -c 2 -o concat > eRNA.f10.VISTA.txt

# ====================================
# DNase signal density 
# ====================================
# fetal brain DNase data from Roadmap (E72, E73)
# https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15

bigWigAverageOverBed $EXTERNAL_FEATURE/DNase/DNase.Roadmap.fBrain.pval.bigwig $inputbed stdout | sort -k1,1 | cut -f1,6 > eRNA.f11.DNase.txt

# ====================================
# DNase cluster 
# ====================================
# download from ENCODE DNase cluster
# ( DNase cluster in V2: minimal score 500, at least 5 cell lines): http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeRegDnaseClustered
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/DNase/wgEncodeRegDnaseClusteredV2.filtered.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f12.DNaseENCODE.txt

# Roadmap DNase (http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation)
#for i in E071 E074 E068 E069 E072 E067 E073 E070 E082 E081; do echo $i; curl -s http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_enh/regions_enh_${i}.bed.gz > externalData/DNase/regions_enh_${i}.bed.gz; done
#zcat externalData/DNase/regions_enh_*.bed.gz | sortBed | mergeBed -i - > externalData/DNase/regions_enh_merged.brain.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/DNase/regions_enh_merged.brain.narrowPeak -c | sort -k4,4 | cut -f4,5 > eRNA.f12.DNaseROADMAP.txt


# ====================================
# Conservation - mean phyloP score
# ====================================
#rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate ../Conservation
#zcat ../Conservation/vertebrate/chr*.gz | gzip -c > ../Conservation/vertebrate/phyloP46way.wigFix.gz
#wigToBigWig ../Conservation/vertebrate/phyloP46way.wigFix.gz $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size ../Conservation/vertebrate/phyloP46way.wigFix.bigwig

bigWigAverageOverBed $EXTERNAL_FEATURE/Conservation/Conservation.phyloP46way.bigwig $inputbed stdout | sort -k1,1 | cut -f1,6 > eRNA.f13.phyloP.txt

# ====================================
# bConserved2zf - conserved to zebrafish
# ====================================
## download liftover chain 
#curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToDanRer7.over.chain.gz | gunzip > externalData/Conservation/hg19ToDanRer7.over.chain
# require a minimum ratio of 10% that must remap
liftOver $inputbed $EXTERNAL_FEATURE/Conservation/hg19ToDanRer7.over.chain stdout unmapped.liftover -minMatch=0.1 | awk '{OFS="\t"; print $4, "zv9|"$1":"$2"-"$3;}' | sort -k1,1 | join -1 4 -2 1 -a 1 -o '0,2.2' -e 'NA' <(sort -k4,4 $inputbed) - | sed 's/ /\t/' > eRNA.f14.bConserved2zf.txt

# ====================================
# Conservation - overlap with HCNE or not
# ====================================
# overlap with HCNE
# curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_hg19_danRer7_70pc_50col.bed.gz | gunzip > ../Conservation/HCNE_hg19_danRer7_70pc_50col.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f15.HCNE.txt

# ====================================
# GWAS (see eRNA.SNP.enrichment.sh for data preparation)
# ====================================
# 1. merge the recent PD-GWAS and NHGRI GWAS catalog SNPs 
# Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.bed (generated by Ganqiang)
# 2. get SNPs in LD
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.PLINK.LD_w250.r2_0.8.bed.autosomal.associations.bed
intersectBed -a $inputbed -b $snps_in_LD -c | sort -k4,4 | cut -f4,5 > eRNA.f16.GWAS.txt
# optional TODO: filter only the brain specific GWAS SNPs, or just PD GWAS

# ====================================
# eQTL of eRNA 
# ====================================

#mkdir eQTL; cd eQTL
#awk 'BEGIN{OFS="\t"; print "locus\tchr\ts1\ts2"}{print $4,$1,$2,$3}' ../$inputbed > eRNA.loci.txt
#ln -fs ../eRNA.90samples.meanRPM.xls eRNA.RPKM.xls
#module unload R/3.1.0; module load R/3.0.2
#bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_PEER_eQTL.eRNA.R
#module unload R/3.0.2; module load R/3.1.0
#cd ..

cut -f2 ../../eQTL/eRNA/final.cis.eQTL.FDR.05.xls | grep chr  | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | join -1 4 -2 2 -a 1 -t $'\t' -e 0 -o '0,2.1' <(sort -k4,4 $inputbed) -  > eRNA.f18.eSNP.txt

# ====================================
# number of HiTNEs in host genes.
# ====================================
# if overlap with multiple genes, take the longest one
intersectBed -a $inputbed -b $GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k12,12nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8,12 | sed 's/\t\./\tNA/g' > $TMPFILE
# host gene (if any); otherwise NA
cut -f4-5 $TMPFILE > eRNA.f19.Hostgene.txt

grep -vw NA $TMPFILE | cut -f5 | sort | uniq -c | join -1 5 -2 2 -a 1 -e '0' -o '1.1,1.2,1.3,1.4,1.5,1.6,2.1' <(sort -k5,5 $TMPFILE) - | sed 's/ /\t/g' > $TMPFILE.2

cut -f4,7 $TMPFILE.2 | sort -k1,1 > eRNA.f20.nHostgene.txt

# ====================================
# length of host genes
# ====================================
cut -f4,6 $TMPFILE.2 | sort -k1,1 > eRNA.f21.lenHostgene.txt

# ====================================
# length of meta intron of host gene (see REAME.txt in $GENOME/Annotation/Genes for how to generate meta intron and exons)
# ====================================
cat $GENOME/Annotation/Genes/exons.meta.bed | awk '{OFS="\t"; print $0,$3-$2}' | sort -k4,4 | groupBy -g 4 -c 7 -o sum > $TMPFILE.3
cat $GENOME/Annotation/Genes/introns.meta.bed | awk '{OFS="\t"; print $0,$3-$2}' | sort -k4,4 | groupBy -g 4 -c 7 -o sum | join -1 1 -2 1 -a 1 -e '0' -o '0,1.2,2.2' $TMPFILE.3 - | sed 's/ /\t/g' > $GENOME/Annotation/Genes/genes.meta.exon.intron.length.txt

join -1 5 -2 1 -a 1 -e '0' -o '1.4,2.3' <(sort -k5,5 $TMPFILE.2) <(sort -k1,1 $GENOME/Annotation/Genes/genes.meta.exon.intron.length.txt) |sed 's/ /\t/g' | sort -k1,1 > eRNA.f22.lenHostgeneMetaintron.txt
 
# ====================================
# overlap with HeLa S3 enhancers (defined by CAGE or histone)
# ====================================
# download the Binary matrix of enhancer usage file, if the enhancer is 1 in HeLa cell line, it should be a HeLa enhancer
# curl -s http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancer_usage.csv.gz > externalData/CAGE/hg19_permissive_enhancer_usage.csv.gz
#zcat externalData/CAGE/hg19_permissive_enhancer_usage.csv.gz | rowsToCols stdin stdout -fs=',' | grep -P "858648|Hela"  | rowsToCols stdin stdout -tab | awk '{OFS="\t"; for(i=1;i<NF;i++) s+=$i; if(s>0) print $NF}' | sed 's/"//g;s/[:-]/\t/g' > externalData/CAGE/hg19_permissive_enhancer.HeLaS3.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/hg19_permissive_enhancer.HeLaS3.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f23.CAGE_HeLa.txt

# Use the enhancers from the combined segmentation (Segway + chromHMM) in HeLaS3 cell line. See reference:
# http://genome-test.cse.ucsc.edu/cgi-bin/hgTrackUi?hgsid=389206706_Fowrtl09k5DnQ80vvZaeHZAP4hkI&c=chr21&g=wgEncodeAwgSegmentation
#curl -s http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedHelas3.bed.gz | gunzip | fgrep E > externalData/Segment/wgEncodeAwgSegmentationCombinedHelas3.Enh.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/Segment/wgEncodeAwgSegmentationCombinedHelas3.Enh.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f24.chromHMM_HeLa.txt

# ====================================
# overlap with SK-N-MC enhancers (defined by CAGE)
# ====================================
# download the Binary matrix of enhancer usage file, if the enhancer is 1 in HeLa cell line, it should be a HeLa enhancer
# curl -s http://enhancer.binf.ku.dk/presets/hg19_permissive_enhancer_usage.csv.gz > externalData/CAGE/hg19_permissive_enhancer_usage.csv.gz
#zcat externalData/CAGE/hg19_permissive_enhancer_usage.csv.gz | rowsToCols stdin stdout -fs=',' | grep -P "858648|SK-N-MC"  | rowsToCols stdin stdout -tab | awk '{OFS="\t"; for(i=1;i<NF;i++) s+=$i; if(s>0) print $NF}' | sed 's/"//g;s/[:-]/\t/g' > externalData/CAGE/hg19_permissive_enhancer.SK-N-MC.bed
intersectBed -a $inputbed -b $EXTERNAL_FEATURE/CAGE/hg19_permissive_enhancer.SK-N-MC.bed -c | sort -k4,4 | cut -f4,5 > eRNA.f25.CAGE_SK-N-MC.txt

# ====================================
# if overlap with any eQTL of genes
# ====================================
fgrep -w -f <(cut -f1 ~/neurogen/rnaseq_PD/results/eQTL/HCILBSNDA89samples/final.cis.eQTL.FDR.05.xls | grep -v SNP | sort -u) /data/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID | awk '{OFS="\t"; print $2,$3-1,$3,$1}' | intersectBed -a $inputbed -b stdin -c | sort -k4,4 | cut -f4,5 > eRNA.f28.eGene.txt


## merge into a big file: eRNA.characterize.xls
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.characterize.merge.R `ls eRNA.f*.txt`

cd ~/eRNAseq

exit;




















# ====================================
### raw reads count
# ====================================
# Note: use '-split' to exclude spliced reads when counting for intronic eRNA
for i in ~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.bam; do bsub -q normal -n 1 $HOME/neurogen/pipeline/RNAseq/modules/_get_readscount_per_region.sh $i $inputbed $i.eRNA.rawcount; done

## read in only the 80 subjects/samples (w/ genotype)
R
path="~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*/uniq/accepted_hits.bam.eRNA.rawcount"
IDs=read.table('~/neurogen/rnaseq_PD/results/merged/RNAseqID.wGenotyped.list',stringsAsFactors =F)[,1]
IDs=IDs[grep("^[HI].*_SNDA", IDs)]
EXP=data.frame(); id="locus"
for(i in Sys.glob(path)){
    ii=sub(".*run_output/(.*)/uniq.*","\\1", i);
    if(! (ii %in% IDs)) next;
    id=c(id, ii);    print(i)
    # read background
    expression=read.table(i, header=F) 
    if(ncol(EXP)==0) EXP=expression[,4:5] else EXP=cbind(EXP, expression[,5]); 
}
colnames(EXP)=id;
write.table(EXP, "eRNA.80samples.rawcount.xls", col.names=T, row.names=F, sep="\t", quote=F)
q('no')



# measure the splicing ratio & intron coverage of its host gene (if any), in order to remove the eRNAs likely being from the fuzzy genes (e.g. highly actively transcribed genes, co-transcriptional splicing? pre-mRNA? if intronic coverage > 30%)
# =================
# run premRNA.sh first
#--------
intersectBed -a $inputbed -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k15,15nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8-13,15 | sed 's/\t\./\t-1/g' > eRNA.premRNAratio

#intersectBed -a $inputbed -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wo | cut -f8 | sort | uniq -c | join -1 5 -2 2 -a 1 -e '0' -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.1' <(sort -k5,5 eRNA.premRNAratio) - | sed 's/ /\t/g' | sortBed > eRNA.premRNAratio.bed
#
#awk '{OFS="\t"; if($5!=-1) print $3-$2,$5,$11,$12;}' eRNA.premRNAratio.bed | sort -k2,2 | groupBy -g 2,3,4 -c 1 -o sum > genes.geneLength.eRNAcount.eRNAlength.tab

# how many eRNA are expected in the null distribution
sortBed -i /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/exons.meta.bed | mergeBed | cat - $GENOME/Annotation/Genes/hg19.gap.bed | bedtools shuffle -excl stdin -i $inputbed -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a - -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wao | awk '{OFS="\t"; print $0,$7-$6;}' | sort -k4,4 -k15,15nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f1-4,8-13,15 | sed 's/\t\./\t-1/g' > eRNA.shuffle.premRNAratio

## R
R
df1=read.table("eRNA.premRNAratio", header=F)
colnames(df1)=c("chr","start","end","ID", "host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length")
df2=read.table("eRNA.shuffle.premRNAratio", header=F)
colnames(df2)=c("chr","start","end","ID", "host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length")

library(ggplot2)
m <- ggplot(df1, aes(x = end-start))
m + geom_density() + scale_x_log10() +
geom_density(data=subset(df1,cov_intron>0.3), aes(x= end-start), color = "blue") +
geom_density(data=subset(df1,cov_intron<=0.3), aes(x= end-start), color = "red")
ggsave("eRNA.size.pdf", width=6, height=6)

require(plyr)
df=ddply(subset(cbind(df1, eRNAlength=df1$end-df1$start), host_gene!="-1"), c("host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length"), summarize, eRNAcount=length(hostGene_length), eRNAlength=sum(eRNAlength))
df=cbind(df, type="HiTNee")
df=rbind(df, cbind(ddply(subset(cbind(df2, eRNAlength=df2$end-df2$start), host_gene!="-1"), c("host_gene","rpkm_exon", "rpkm_intron", "cov_exon", "cov_intron", "splicing_ratio", "hostGene_length"), summarize, eRNAcount=length(hostGene_length), eRNAlength=sum(eRNAlength)), type="shuffled"))

library(ggplot2)
ggplot(df, aes(eRNAcount, fill = type)) + geom_density(alpha = 0.2,adjust=3) + scale_x_log10() + theme_bw()+ theme(legend.position=c(1, 1), legend.justification=c(1,1))
ggsave("eRNAcount.in.hostgene.pdf", width=6, height=4)

df0=subset(df, host_gene %in% subset(df, type=='HiTNee' & eRNAcount>100)$host_gene)
df0$host_gene=gsub()
df0$host_gene2 = reorder(df0$host_gene, -df0$eRNAcount)
ggplot(df0, aes(x=host_gene2, y=eRNAcount, fill=type)) +
geom_bar(stat = "identity", position="dodge") +
coord_flip()+
theme_bw()+
theme(legend.position=c(1, 1), legend.justification=c(1,1), axis.text.y = element_text(size=5))

barplot(eRNAcount~type, df=subset(df, host_gene %in% subset(df, type=='HiTNee' & eRNAcount>100)$host_gene))
, aes(eRNAcount, fill=type)) + geom_bar() + theme_bw()+ theme(legend.position=c(1, 1), legend.justification=c(1,1))

df=df[with(df,order(type, -eRNAcount)),]
df=read.table("genes.eRNAcount.eRNAlength.tab")
colnames(df)=c("geneID","eRNAcount","geneLength")

pdf("eRNA.hostgene.pdf", width=12, height=8)
plot(eRNAcount~geneLength, df, col=ifelse(grepl("protein_coding",geneID),'#00000066','#ff000066'), pch=20)
text(df$geneLength, df$eRNAcount, labels=ifelse(df$eRNAcount>100,gsub("___ENSG.*","",df$geneID),""), cex=0.5, pos=4, offset=0.2, col=ifelse(grepl("protein_coding",df$geneID),'#000000','#ff0000'))

df=subset(df, eRNAcount>100)
df=df[with(df, order(-eRNAcount)),]
barplot(df$eRNAcount,
    names.arg =gsub("___ENSG.*","",df$geneID),
    ylab="eRNAs count in host gene",
    las=2, cex.names=.5,
    col=ifelse(grepl("protein_coding",df$geneID),"gray","red"),
    space=0.8,
    legend.text = c("protein-coding", "non-protein-coding"),
    args.legend = list(x = "topright", bty = 'n', fill=c("gray","red")))

dev.off()

q('no')




# fisher exact test to see the signifiance

R
df=read.table("eRNA.overlap.txt")
rownames(df)=df[,1]; df=df[,-1]
barplot(table(rowSums(df>0)), log='y', axes=F)
axis(2, 0:5, labels=10**c(0:5))

## robust eRNAs: at least one line of evidences from CAGE|Histone|TFBS|HCNE|eQTL|GWAS
echo -e "ID\tHistone\tCAGE\tDNase\tTFBS\tHCNE\te2QTL\teQTL\tGWAS\tn_evidences\tRPKM\thost_gene\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > robust_eRNAs.info.tab
awk '{s=0; for(i=2;i<=NF;i++) if($i>0) s++; OFS="\t"; print $0,s}' eRNA.overlap.txt | paste - <(cut -f2 eRNA.RPKM) | paste - <(cut -f5- eRNA.premRNAratio) | awk '$10>=3 && ($2>0 || $3>0 || $5>0 || $6>0 || $8>0 || $9>0)' >> robust_eRNAs.info.tab

echo -e "ID\tHistone\tCAGE\tDNase\tTFBS\tHCNE\te2QTL\teQTL\tGWAS\tn_evidences\tRPKM\thost_gene\trpkm_exon\trpkm_intron\tcov_exon\tcov_intron\tsplicing_ratio" > permissive_eRNAs.info.tab
awk '{s=0; for(i=2;i<=NF;i++) if($i>0) s++; OFS="\t"; print $0,s}' eRNA.overlap.txt | paste - <(cut -f2 eRNA.RPKM) | paste - <(cut -f5- eRNA.premRNAratio) | awk '$10>0' >> permissive_eRNAs.info.tab

## get UCSC screenshot for eRNAs
grep chr robust_eRNAs.info.tab | cut -f1 | awk '{OFS="\t"; split($1,a,"_");print a[1],a[2],a[3],$1}' | Rscript ~/neurogen/pipeline/RNAseq/modules/_UCSC2PDF.R stdin
Rscript ~/neurogen/pipeline/RNAseq/modules/_table2html.R robust_eRNAs.info.tab robust.html robust_eRNAs


grep chr permissive_eRNAs.info.tab | cut -f1,10 | awk '{OFS="\t"; split($1,a,"_");print a[1],a[2],a[3],$1,$2}' | intersectBed -a stdin -b gencode.v19.longestTx.exons-introns.RNAseq.aboveBasal.bigwig.bed -wo | awk '{OFS="\t"; print $0,$8-$7;}' |sort -k4,4 -k16,16nr | awk '{OFS="\t"; if($4!=id) {print; id=$4;}}' | cut -f4,5,9,16 > eRNA.hostgene.tab

R
df=read.table("eRNA.hostgene.tab")
colnames(df)=c('ID','n_evidence','hostGene','length')
library("plyr")
df2=ddply(df, .(hostGene, length), summarize, n_eRNA=length(ID),sum_evidence=sum(n_evidence>=3))
df=ddply(df2, .(n_eRNA), summarize, n_gene=length(hostGene), mean_length=mean(length), mean_evidence=mean(sum_evidence))
pdf("eRNA.hostgene.pdf", width=10, height=4)
plot(df$n_gene, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
plot(df$mean_length, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
plot(df$mean_evidence, type='h', xaxt='n'); axis(1, at=1:nrow(df), labels=df$n_eRNA, las=2, cex.axis=0.7)
dev.off()



#############################################################
# measure eRNA, binarily and continously, with other features (e.g. CAGE, DNase, histone etc.)
#############################################################
# CAGE+
~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../CAGE/CAGE.fwd.bigwig $inputbed 1 max | sort -k1,1 > eRNA.onOtherFeatures.txt

# CAGE-
~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../CAGE/CAGE.rev.bigwig $inputbed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
mv tmp.list eRNA.onOtherFeatures.txt

# H3k4me1, H3K4me3 H3K27ac K3K27me3 H3K36me3 H3K9ac
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K9ac; do
    echo $i;
    ~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../Histone/Histone.SN.$i.bigwig $inputbed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
    mv tmp.list eRNA.onOtherFeatures.txt
done

# DNase, TFBS, Conservation
for i in DNase TFBS Conservation; do
    echo $i;
    ~/neurogen/pipeline/RNAseq/bin/toBinRegionsOnBigwig.sh ../$i/$i.bigwig $inputbed 1 max | sort -k1,1 | cut -f2 | paste eRNA.onOtherFeatures.txt - > tmp.list
    mv tmp.list eRNA.onOtherFeatures.txt
done

# add header
echo -e "locus\tCAGE.fwd\tCAGE.rev\tH3K4me1\tH3K4me3\tH3K27ac\tH3K27me3\tH3K36me3\tH3K9ac\tDNase\tTFBS\tConservation" > tmp.list
cat eRNA.onOtherFeatures.txt >> tmp.list
mv tmp.list eRNA.onOtherFeatures.txt

## clustering
R
df=read.table("eRNA.onOtherFeatures.txt", header=T)
rownames(df)=df[,1]; df=df[,-1]

cage=df[,c('CAGE.fwd','CAGE.rev')]
library(pheatmap) # IF NOT, install.packages('pheatmap')
library("RColorBrewer")

pheatmap(log(1+df),scale='column',
    show_rownames=F,
    cluster_cols=F,
    clustering_distance_rows="correlation",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    filename="eRNA.onOtherFeatures.png")
)
