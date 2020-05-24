## Similar as eRNA.SNP.enrichment.sh, but rather use all SNPs in LD, just use GWAS and PICS SNPs (similar as Fig 2a in PMID: 26780995)
# Usage: 
# for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 
# for i in HCILB_SNDA HC_TCPY HC_MCPY HC_FB HC_PBMC; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 

# alternatively,
## for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh PLINK $i; done 
# bsub -n 1 -q normal -J HCILB_SNDA bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP HCILB_SNDA

type=$1
samplegroup=$2

# debug
# type='GWAS-DA'; samplegroup='HCILB_SNDA'

cd ~/eRNAseq/$samplegroup

## pre-steps to get GWAS SNPs and PICS SNPs: 
## Ref ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/README.txt

## 8 ND diseases: AD, FTD, MSA, myasthenia gravis (MG), Charcot Marie Tooth (CMT), PSP, ALS, and PD. (from NeuroX paper)
# awk '{FS="\t";OFS="\t"; print $7,$0;}' gwas_catalog_v1.0-downloaded.hg19.bed  | grep -iE "^Alzheimer's disease|^amyotrophic lateral sclerosis|^multiple sclerosis|^Parkinson|^frontotemporal dementia|^multiple system atrophy|^myasthenia gravis|^Charcot Marie Tooth|^progressive supranuclear palsy" | cut -f2-8 > gwas_catalog_v1.0-downloaded.hg19.ND.bed
## 11 DA (dopamine) diseases/traits: response to antipsychotic treatment|sleep quality, insomnia, sleep-related phenotypes, schizophrenia, PD, ADHD, response to methylphenidate treatment, addiction, bipolar disorder, Response to iloperidone treatment
# awk '{FS="\t";OFS="\t"; print $7,$0;}' gwas_catalog_v1.0-downloaded.hg19.bed  | grep -iE "response to antipsychotic treatment|sleep quality|insomnia|sleep-related phenotypes|schizophrenia|Parkinson|attention deficit hyperactivity disorder|response to methylphenidate treatment|addiction|bipolar disorder|^Response to iloperidone treatment" | cut -f2-8 > gwas_catalog_v1.0-downloaded.hg19.DA.bed

[ "$type" == "GWAS" ] && snps=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pruned.bed 
[ "$type" == "GWAS-ND" ] && snps=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.ND.bed
[ "$type" == "GWAS-DA" ] && snps=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.DA.bed
[ "$type" == "PICS" ] && snps=$GENOME/Annotation/GWASCatalog/PICS_downloaded_2017-12-11.hg19.bed

## extract all autosomal.associations
awk 'BEGIN{FS="\t"; OFS="\t";}{print $1,$2,$3,$7;}' $snps | grep -v chrX | grep -v chrY | sort -u > $snps.autosomal.associations.bed

#### prepare data
## random regions
#bedtools shuffle -excl <(cat ../blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed) -noOverlapping -i eRNA.bed -g $GENOME/Annotation/Genes/ChromInfo.txt > eRNA.random.Lmatched.bed
#bedtools random -n 100000 -l 400 -seed 1234 -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(cat ../blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -v > eRNA.random.L400.bed
ln -fs eRNA.random.L400.bed eRNA.random.bed

## mRNA inner exons
grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed | awk '{if(id!=$4) id=$4; else print}' | sort -k4,4 -k2,2nr | awk '{if(id!=$4) id=$4; else print}' > $GENOME/Annotation/Genes/mRNA.innner.exon.bed

## promoter -- [-200,+200] of protein-coding GENCODE v19 TSS
grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400}' > $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed

### overlapping 

EXTERNAL_FEATURE=~/eRNAseq/externalData
# total
echo total `wc -l $snps | cut -f1 -d' '` `grep -Ev "_|M|X|Y" $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size | awk '{s+=$2}END{print s}'` > SNP.$type.counts.summary
# 5UTR
echo utr5 `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/utr5.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/utr5.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# innner exon
echo coding_exon `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/mRNA.innner.exon.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/mRNA.innner.exon.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# 3UTR
echo utr3 `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/utr3.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/utr3.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# intron
echo intron `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/introns.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/introns.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# intergenic
echo intergenic `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/intergenic.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/intergenic.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# promoter
echo promoter `intersectBed -a $snps.autosomal.associations.bed -b $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed -u | wc -l` `sortBed -i $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# FANTOM5-enhancer
echo FANTOM5_enhancer `intersectBed -a $snps.autosomal.associations.bed -b $EXTERNAL_FEATURE/CAGE/permissive_enhancers.bed -u | wc -l` `awk '{s+=($3-$2);}END{print s}' $EXTERNAL_FEATURE/CAGE/permissive_enhancers.bed` >> SNP.$type.counts.summary
# chromHMM-enhancer-cellline
echo chromHMM_enhancer `intersectBed -a $snps.autosomal.associations.bed -b $EXTERNAL_FEATURE/Segment/wgEncodeBroadHmm.strongEnhancer.bed -u | wc -l` `sortBed -i $EXTERNAL_FEATURE/Segment/wgEncodeBroadHmm.strongEnhancer.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# chromHMM-enhancer-brain
echo chromHMM_enhancer_brain `intersectBed -a $snps.autosomal.associations.bed -b $EXTERNAL_FEATURE/Segment/15_coreMarks_segments.E6E7E12.bed -u | wc -l` `sortBed -i $EXTERNAL_FEATURE/Segment/15_coreMarks_segments.E6E7E12.bed | mergeBed | awk '{s+=($3-$2);}END{print s}'` >> SNP.$type.counts.summary
# TNE
echo TNE `intersectBed -a $snps.autosomal.associations.bed -b eRNA.bed -u | wc -l` `awk '{s+=($3-$2);}END{print s}' eRNA.bed` >> SNP.$type.counts.summary
# TNE-random
echo random `intersectBed -a $snps.autosomal.associations.bed -b eRNA.random.bed -u | wc -l` `awk '{s+=($3-$2);}END{print s}' eRNA.random.bed` >> SNP.$type.counts.summary

#echo "## Make bar plot" 
### ##################
Rscript $pipeline_path/src/eRNA.SNP.enrichment.PICS.R $type $samplegroup