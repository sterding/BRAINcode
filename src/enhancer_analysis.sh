## script to analyse eRNA defined by HCILB RNA-seq

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

cd ~/projects/PD/results/eRNA/externalData/RNAseq

################################################################################################
# define eRNA
################################################################################################
for i in HCILB_SNDA HC_nonNeuron HC_PY; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.define.sh $i; done
for i in HC_FB HC_PBMC HC_MCPY HC_TCPY; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.define.sh $i; done
for i in ILB_SNDA PD_SNDA; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.define.sh $i; done

# use FDR as correction method for all cell types
for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY ILB_SNDA PD_SNDA; do cd $i; ln -fs eRNA.fdr.bed eRNA.bed; cd -; done

## make available to UCSC
cd ~/neurogen/rnaseq_PD/for_display/; bash $pipeline_path/src/_make_trackDb.HTNE.sh ; cd -

################################################################################################
# define shared and private eRNA
################################################################################################
# ven diagram for overlap of 3 major groups (a: SNDA; b:PY; c:NN)
echo a `intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_PY/eRNA.bed HC_nonNeuron/eRNA.bed) -v | wc -l` > venn.diagram.fdr.txt
echo b `intersectBed -a HC_PY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_nonNeuron/eRNA.bed) -v | wc -l` >> venn.diagram.fdr.txt
echo c `intersectBed -a HC_nonNeuron/eRNA.bed -b <(cat HC_PY/eRNA.bed HCILB_SNDA/eRNA.bed) -v | wc -l` >> venn.diagram.fdr.txt
echo ab `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_PY/eRNA.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.bed -v | wc -l` >> venn.diagram.fdr.txt
echo ac `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -v | wc -l` >> venn.diagram.fdr.txt
echo bc `intersectBed -a HC_PY/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bed -v | wc -l` >> venn.diagram.fdr.txt
echo abc `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -u | wc -l` >> venn.diagram.fdr.txt

echo a `intersectBed -a HCILB_SNDA/eRNA.bonferroni.bed -b <(cat HC_PY/eRNA.bonferroni.bed HC_nonNeuron/eRNA.bonferroni.bed) -v | wc -l` > venn.diagram.bonferroni.txt
echo b `intersectBed -a HC_PY/eRNA.bonferroni.bed -b <(cat HCILB_SNDA/eRNA.bonferroni.bed HC_nonNeuron/eRNA.bonferroni.bed) -v | wc -l` >> venn.diagram.bonferroni.txt
echo c `intersectBed -a HC_nonNeuron/eRNA.bonferroni.bed -b <(cat HC_PY/eRNA.bonferroni.bed HCILB_SNDA/eRNA.bonferroni.bed) -v | wc -l` >> venn.diagram.bonferroni.txt
echo ab `intersectBed -a HCILB_SNDA/eRNA.bonferroni.bed -b HC_PY/eRNA.bonferroni.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.bonferroni.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo ac `intersectBed -a HCILB_SNDA/eRNA.bonferroni.bed -b HC_nonNeuron/eRNA.bonferroni.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bonferroni.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo bc `intersectBed -a HC_PY/eRNA.bonferroni.bed -b HC_nonNeuron/eRNA.bonferroni.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bonferroni.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo abc `intersectBed -a HCILB_SNDA/eRNA.bonferroni.bed -b HC_nonNeuron/eRNA.bonferroni.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bonferroni.bed -u | wc -l` >> venn.diagram.bonferroni.txt

# celltype-private ones
intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_PY/eRNA.bed HC_nonNeuron/eRNA.bed) -v > HCILB_SNDA/eRNA.private.major.bed
intersectBed -a HC_PY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_nonNeuron/eRNA.bed) -v > HC_PY/eRNA.private.major.bed
intersectBed -a HC_nonNeuron/eRNA.bed -b <(cat HC_PY/eRNA.bed HCILB_SNDA/eRNA.bed) -v > HC_nonNeuron/eRNA.private.major.bed

intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_PY/eRNA.bed HC_nonNeuron/eRNA.bed) -v | awk '{OFS="\t"; print $0,"SNDAonly"}' > HCILB_SNDA/eRNA.privacy.bed
intersectBed -a HCILB_SNDA/eRNA.bed -b HC_PY/eRNA.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.bed -v | awk '{OFS="\t"; print $0,"SNDA|PY"}' >> HCILB_SNDA/eRNA.privacy.bed
intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -v | awk '{OFS="\t"; print $0,"SNDA|NN"}' >> HCILB_SNDA/eRNA.privacy.bed
intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -u | awk '{OFS="\t"; print $0,"SNDA|PY|NN"}' >> HCILB_SNDA/eRNA.privacy.bed

intersectBed -a HC_PY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_nonNeuron/eRNA.bed) -v | awk '{OFS="\t"; print $0,"PYonly"}' > HC_PY/eRNA.privacy.bed
intersectBed -a HC_PY/eRNA.bed -b HCILB_SNDA/eRNA.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.bed -v | awk '{OFS="\t"; print $0,"PY|SNDA"}' >> HC_PY/eRNA.privacy.bed
intersectBed -a HC_PY/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bed -v | awk '{OFS="\t"; print $0,"PY|NN"}' >> HC_PY/eRNA.privacy.bed
intersectBed -a HC_PY/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bed -u | awk '{OFS="\t"; print $0,"SNDA|PY|NN"}' >> HC_PY/eRNA.privacy.bed

intersectBed -a HC_nonNeuron/eRNA.bed -b <(cat HC_PY/eRNA.bed HCILB_SNDA/eRNA.bed) -v | awk '{OFS="\t"; print $0,"NNonly"}' > HC_nonNeuron/eRNA.privacy.bed
intersectBed -a HC_nonNeuron/eRNA.bed -b HC_PY/eRNA.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bed -v | awk '{OFS="\t"; print $0,"NN|PY"}' >> HC_nonNeuron/eRNA.privacy.bed
intersectBed -a HC_nonNeuron/eRNA.bed -b HCILB_SNDA/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -v | awk '{OFS="\t"; print $0,"NN|SNDA"}' >> HC_nonNeuron/eRNA.privacy.bed
intersectBed -a HC_nonNeuron/eRNA.bed -b HCILB_SNDA/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -u | awk '{OFS="\t"; print $0,"SNDA|PY|NN"}' >> HC_nonNeuron/eRNA.privacy.bed

# for minor cell tyes
intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HC_MCPY/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HCILB_SNDA/eRNA.private.minor.bed
intersectBed -a HC_TCPY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_MCPY/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_TCPY/eRNA.private.minor.bed
intersectBed -a HC_MCPY/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_MCPY/eRNA.private.minor.bed
intersectBed -a HC_FB/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_MCPY/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_FB/eRNA.private.minor.bed
intersectBed -a HC_PBMC/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_FB/eRNA.bed HC_MCPY/eRNA.bed) -v > HC_PBMC/eRNA.private.minor.bed

# shared btw SNDA and the other
> HCILB_SNDA/eRNA.shared.minor.txt
for i in HC_TCPY HC_MCPY HC_PBMC HC_FB; do echo $i; echo "HCILB_SNDA vs. $i:" `intersectBed -a HCILB_SNDA/eRNA.bed -b $i/eRNA.bed -u | wc -l` >> HCILB_SNDA/eRNA.shared.minor.txt; done

## fold change and p-value for private HTNEs
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HCILB_SNDA/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HC_PY/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HC_nonNeuron/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA

cut -f4 eRNA.private.major.bed | fgrep -f - eRNA.meanRPM.allSamples.ttest.txt | awk '$2>=2 && $3<0.01 && $4>=2 && $5<0.01' | cut -f1 | fgrep -f - eRNA.characterize.xls | awk '$28==1 && $6>=10 && $7>0 && $8>0' 
chr1_66277655_66278270	19766	90.8832	0.133147	0.0339962	36	1	1	1	2	NA	0.622140.185885	NA	0	0	0	PDE4B___ENSG00000184588.13___protein_coding	223	582063	573422	1	1	1
chr18_13566745_13567739	349746	93.3733	0.170893	0.305512	50	2	2	2	10	NA	1.23605-0.0844749	NA	0	0	0	LDLRAD4___ENSG00000168675.14___protein_coding	144	435258	418499	2	1
chr8_41570144_41570663	183877	65.0584	0.0954287	0.505082	32	2	1	0	9	NA	40.64260.390697	NA	0	0	0	ANK1___ENSG00000029534.15___protein_coding	63	243542	232743	0	0	1
cat eRNA.meanRPM.allSamples.ttest.txt | awk '$2>=2 && $3<0.01 && $4>=2 && $5<0.01' | cut -f1 | fgrep -f - eRNA.characterize.xls | awk '$28==1 && $6>=10 && $7>0 && $8>0'  # n=4
cat eRNA.meanRPM.allSamples.ttest.txt | awk '$2>=1 && $3<0.01 && $4>=1 && $5<0.01' | wc -l  # n=10295
cut -f4 eRNA.private.major.bed | fgrep -f - eRNA.meanRPM.allSamples.ttest.txt | awk '$2>=1 && $3<0.01 && $4>=1 && $5<0.01' | wc -l  # n=9777 (95%)



join -1 4 -2 1 <(sort -k4,4 eRNA.privacy.bed) <(sort -k1,1 eRNA.meanRPM.allSamples.ttest.txt) > eRNA.privacy.ttest.txt
R
setwd("~/eRNAseq/HCILB_SNDA")
df=read.table("eRNA.privacy.ttest.txt", header=F, check.names = F)
head(df); str(df)

pdf("eRNA.privacy.ttest.pdf", width=5, height=5)
plot(df$V6, -log10(df$V7), bg=adjustcolor(as.integer(df$V5), 0.8), pch=21,cex=.7, col='white',xlab='log2 fold-change', ylab='-log10(p-value)', main='SNDA vs. PY')
legend("topleft", legend=levels(df$V5), col=adjustcolor(1:length(levels(df$V5)), 0.8), pch=20, cex=1)

plot(df$V8, -log10(df$V9), bg=adjustcolor(as.integer(df$V5), 0.8), pch=21,cex=.7, col='white',xlab='log2 fold-change', ylab='-log10(p-value)', main='SNDA vs. NN')
legend("topleft", legend=levels(df$V5), col=adjustcolor(1:length(levels(df$V5)), 0.8), pch=20, cex=1)

dev.off()

#############################################################
# GWAS SNP enrichment
#############################################################

for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh PLINK $i; done 
for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 
Rscript $pipeline_path/src/eRNA.SNP.enrichment.merge.R HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY

# private only
echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.private.major.counts.summary
echo "HCILB_SNDA" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b HCILB_SNDA/eRNA.private.bed -u | wc -l | cut -f1 -d' '` >> SNP.private.major.counts.summary
echo "HC_PY" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b HC_PY/eRNA.private.bed -u | wc -l | cut -f1 -d' '` >> SNP.private.major.counts.summary
echo "HC_nonNeuron" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b HC_nonNeuron/eRNA.private.bed -u | wc -l | cut -f1 -d' '` >> SNP.private.major.counts.summary

echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.private.minor.counts.summary
echo "HCILB_SNDA" `cat HCILB_SNDA/eRNA.private.minor.bed | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.private.minor.counts.summary
echo "HC_TCPY" `cat HC_TCPY/eRNA.private.minor.bed | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.private.minor.counts.summary
echo "HC_MCPY" `cat HC_MCPY/eRNA.private.minor.bed | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.private.minor.counts.summary
echo "HC_FB" `cat HC_FB/eRNA.private.minor.bed | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.private.minor.counts.summary
echo "HC_PBMC" `cat HC_PBMC/eRNA.private.minor.bed | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b - -u | wc -l | cut -f1 -d' '` >> SNP.private.minor.counts.summary

Rscript $pipeline_path/src/eRNA.SNP.enrichment.privateOnly.R SNAP minor

## overlap of GWAS --> eQTL --> TFBS 
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
grep Parkinson $snps_in_LD.autosomal.associations.bed | intersectBed -b HCILB_SNDA/eRNA.bed -a - -wo | cut -f1-4,8 | intersectBed -a - -b <(awk '{OFS="\t"; if(NR>1) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2,$5,$6}}' HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.xls) -wao |  sed 's/ /__/g' | sort -k1,1 -k2,2 -k11,11g | awk '{if(s!=$2){print; s=$2}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wao


################################################################################################
# eQTL of eRNA
################################################################################################
cd HCILB_SNDA;
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.R # for HCILB_SNDA only
## run permutations in bash
for i in `seq 1 10000`; do echo $i; [ -e permutations/permutation$i.txt ] || bsub -n 1 -M 500 -q short -J $i Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_permutation_minP.R $i data.RData expression.postSVA.xls 1000; done
# post-eQTL analysis
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.R
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.plot.R rs1684902:61633127:A:G_A:G chr10_61632545_61633312

cat ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.xls | awk '$6<=0.05' > ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.new.d1e6.p1e-2.FDRpt5.xls

## eQTL analysis
# HTNE w/ eQTL
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print;}' final.cis.eQTL.new.d1e6.p1e-2.xls | cut -f2 | sort | uniq -c | sort -k2,2 |awk '{OFS="\t"; print $2,$1}' | join -a 1 -1 1 -2 1 -o "0,1.2,2.6,2.20,2.28" - <(sort eRNA.characterize.xls) | sed 's/\s/\t/g' | cut -f5 | sort | uniq -c
# HTNE w/o eQTL
fgrep -v -f <(awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print;}' final.cis.eQTL.new.d1e6.p1e-2.xls | cut -f2 | sort -u) eRNA.characterize.xls | cut -f28 | sort | uniq -c
# Fisher' test in R
fisher.test(matrix(c(39,23586,106,47291),nrow=2,byrow=T),'greater')  # p-value = 0.1123

## GO analysis for the genes harboring eQTL HTNEs
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print;}' final.cis.eQTL.new.d1e6.p1e-2.xls | cut -f2 | sort | uniq -c | sort -k2,2 |awk '{OFS="\t"; print $2,$1}' | join -a 1 -1 1 -2 1 -o "0,1.2,2.6,2.20,2.28" - <(sort eRNA.characterize.xls) | sed 's/\s/\t/g' | cut -f4 | sort -u | sed 's/___/\t/g' | cut -f2 | sed 's/\..*//g'

## eQTL filters:
## A: colocalize --> FDR --> TFBS/GWAS
## ----------------------------------
awk '{split($2,a,"_"); if(NR==1 || ($7>(a[2]-500) && $7<(a[3]+500))) print}' final.cis.eQTL.new.d1e6.p1.xls  > final.cis.eQTL.new.d500.p1.xls
## final.cis.eQTL.new.d500.p1.xls --> final.cis.eQTL.new.d500.adjusted.xl in R
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' final.cis.eQTL.new.d500.adjusted.xls | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo ## one
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' final.cis.eQTL.new.d500.adjusted.xls | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo ## NONE

## B: GWAS --> FDR --> TFBS
## ----------------------------------
#1. GWAS disease associated
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sort -u > $snps_in_LD.autosomal.associations.bed
more +2 final.cis.eQTL.xls | awk '{OFS="\t"; split($2,a,"_"); print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | cut -f1-8,12 > final.cis.eQTL.GWAS.xls
## out of 495085 GWAS associated SNPs, 440565 (89%) are covered by our SNP sets (genotyping + imputation)

# 2. significant after multi-test correction (FDR)
R
setwd("~/eRNAseq/HCILB_SNDA")
df=read.table("final.cis.eQTL.GWAS.xls", header=F, stringsAsFactors =F, sep="\t", quote="\"")
colnames(df) = c('chr','start','end','SNP','gene','beta','t.stat','p.value','trait')
df$FDR=p.adjust(df$p.value,method='fdr'); df$bonferroni=p.adjust(df$p.value,method='bonferroni')
write.table(df, "final.cis.eQTL.GWAS.adjusted.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)
df=subset(df, FDR<=0.05) # filtered by FDR 
#df=read.table("final.cis.eQTL.GWAS.FDR5pt.xls", header=F, stringsAsFactors =F, sep="\t", quote="\"")
#colnames(df) = c('chr','start','end','SNP','gene','beta','t.stat','p.value','trait','FDR','bonferroni')

pdf("eRNA.eQTL.enrichment.pdf", width=6, height = 8)
# to plot the number of eQTL SNPs per disease
n_eQTL = sort(table(df$trait), decreasing =T)
n_eQTL = n_eQTL[n_eQTL>=3]
length(n_eQTL)
# remove parts in "(*)"
names(n_eQTL) = gsub(" \\(.*","", names(n_eQTL))
par(mar=c(30,5,2,2));
barplot(n_eQTL, width=.2, space=1, las=2, cex=.6,cex.axis=.8, cex.names = .5, col='#F22A7B', border='#F22A7B', ylab="Number of disease eQTL SNPs")
# Manhantan plot for the eQTL, color by diseases
head(df)
rcMpPvals=subset(df, select=c(chr, end, p.value, trait))
chrMax=0;
for (i in 1:22){ 
  ndx <- which(rcMpPvals[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  rcMpPvals[ndx, 2] <- rcMpPvals[ndx, 2] + chrMax;
  chrMax = max(rcMpPvals[ndx, 2]);
  #print(c(i, length(ndx), chrMax))
}
bpMidVec <-c(); labels=c();
for (i in 1:22){
  ndx <- which(rcMpPvals[, 1]==paste0('chr',i))
  if(length(ndx)==0) next
  posSub <- rcMpPvals[ndx, 2]
  bpMidVec <- c(bpMidVec, ((max(posSub) - min(posSub))/2) + min(posSub))
  labels = c(labels, paste0('chr',i))
}

library('ggplot2')
ggplot(rcMpPvals) +
  geom_point(aes(x=end, y=-log10(p.value), colour=as.factor(chr)),size=3, alpha=1/2) +
  #scale_color_manual(values=rep(c('black', 'dark green'), 6)) 
  theme_bw(base_size=15) + 
  theme(legend.position='none') +
  scale_x_continuous(labels=labels, breaks=bpMidVec) +
  #geom_hline(y=4.08, linetype=1, col='red', lwd=1.5) +
  ggtitle('eQTL with disease-associated SNPs') + xlab('') + ylab('-log10(P)')
print(p)
# enrichment test

dev.off()
quit('no')

more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$10<=0.05' > final.cis.eQTL.GWAS.FDR5pt.xls
#more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$11<=0.05' # only have 36 eQTL remained after Bonferroni correction, choose to use FDR

## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.xls # N=1

# 3. disrupt the TFBS motif
TFBS=../externalData/TFBS/factorbookMotifPos.v3.bed
intersectBed -a final.cis.eQTL.GWAS.FDR5pt.xls -b $TFBS -wo > final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls # N=0

cat final.cis.eQTL.GWAS.FDR5pt.TFBS.xls | cut -f9 | sort | uniq -c | sort -k1,1n
grep Parkinson final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

## C: FDR > colocalize --> TFBS: N=0
## ----------------------------------
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05 && $7>=a[2] && $7<=a[3]) print a[1],$7-1,$7,$1,$2;}' final.cis.eQTL.new.d1e6.p1e-2.xls  | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo  # NONE

## D: p<0.01 --> colocalize --> TFBS: N=5
## ----------------------------------
awk '{OFS="\t"; split($2,a,"_"); if($7>=a[2] && $7<=a[3]) print;}' final.cis.eQTL.new.d1e6.p1e-2.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -u  # N=6
awk '{OFS="\t"; split($2,a,"_"); if($7>=a[2] && $7<=a[3]) print;}' final.cis.eQTL.new.d1e6.p1e-2.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<2 || $9>=5) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo > final.cis.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.xls

more +2 final.cis.eQTL.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05 && $7>=a[2] && $7<=a[3]) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo
more +2 final.cis.eQTL.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<3 || $9>10) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo

## disease with significant cooccurance of HTNE-eQTL and GWAS 
## ========================================================
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
# prune GWAS SNPs
more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.txt
more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.txt

R
setwd("~/eRNAseq/HCILB_SNDA")
N=as.numeric(system("wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '", intern=T)) ## total SNPs in dbSNP
#N=as.numeric(system("cut -f1-3 $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed.autosomal.associations.bed | sort -u | wc -l | cut -f1 -d' '", intern=T)) ## total dSNPs in gwas+proxy
n=as.numeric(system("more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '$6<=0.05' | cut -f1,7 | sort -u | wc -l", intern=T))  ## total eSNP
df1=read.table("dSNP.with.eSNP.txt", header=F); rownames(df1)=df1[,1]
df2=read.table("dSNP.without.eSNP.txt", header=F); rownames(df2)=df2[,1]
df=cbind(df1, df2[rownames(df1),2]); df=df[,-1]; colnames(df)=c('dSNP_eSNP','dSNP_N_eSNP')  ## only the disease with both eSNP and dSNP
results = cbind(Disease_or_Trait=rownames(df), 
                df, 
                pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$p.value), 
                OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$estimate))
results = results[with(results, order(-dSNP_eSNP)), ]
results$Disease_or_Trait=gsub("_"," ", results$Disease_or_Trait)
write.table(results, "dSNP.eQTL.enrichment.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)

pdf("dSNP.eQTL.enrichment.pdf", width=3, height=6)
results= subset(results, OR>1 & pvalue<0.01/nrow(df) & dSNP_eSNP>3)  # Bonferroni correction for number of diseases that have both eSNP and dSNP
results$pvalue[results$pvalue==0]=2.2e-16
results$Disease_or_Trait <- factor(results$Disease_or_Trait, unique(as.character(results$Disease_or_Trait)))
library('ggplot2')
p = ggplot(results, aes(x=Disease_or_Trait, y=dSNP_eSNP, fill='#F22A7B',ymax=max(dSNP_eSNP)*1.3)) 
p = p + geom_bar(width=.2, position = position_dodge(width=1), stat="identity") 
p = p + theme_bw() + ylab("Number of eSNP and dSNP") + scale_y_log10()
p = p + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), legend.position='none') 
p = p + geom_text(aes(label=paste0(" (",ceiling(-log10(pvalue)),")")), position = position_dodge(width=1), hjust=0, vjust=.5, angle = 90, size=2.5) 
p = p + ggtitle(paste0(basename(getwd()), " -- dSNP enrichments")) 
print(p)
dev.off()
quit('no')

################################################################################################
# characterize eRNA
################################################################################################
# see eRNA.characterize.sh
for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY ILB_SNDA PD_SNDA; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.characterize.sh $i; done

for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY ILB_SNDA PD_SNDA; do cd $i; Rscript $pipeline_path/src/eRNA.characterize.merge.R `ls eRNA.f*.txt`; cd -; done


cd HCILB_SNDA

# if the overlap is significant or not?
for i in `seq 1 1000`;
do
    echo $i;
    bsub -J random_overlap -oo _random_overlap.$i.log -eo _random_overlap.$i.log -q normal -n 1 -M 500 $pipeline_path/src/overlap.with.random.sh
done

tmp=eRNA.bed
# TFBS hotspot
awk '$6>=1' ~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -c | awk '$5>=5' | wc -l > eRNA.overlap.txt
# P300
awk '$4=="EP300"' ~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -u | wc -l >> eRNA.overlap.txt
# FANTOMF CAGE
intersectBed -a $tmp -b ~/eRNAseq/externalData/CAGE/permissive_enhancers.bed -u | wc -l >> eRNA.overlap.txt
# Roadmpa Histone
intersectBed -a $tmp -b ~/eRNAseq/externalData/Segment/15_coreMarks_segments.E6E7E12.bed -u | wc -l >> eRNA.overlap.txt
# VISTA positive
grep -v NULL ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | intersectBed -a $tmp -b stdin -u | wc -l >> eRNA.overlap.txt
# HCNE
intersectBed -a $tmp -b ~/eRNAseq/externalData/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | wc -l >> eRNA.overlap.txt
# DNase peak from Roadmap
intersectBed -a $tmp -b ~/eRNAseq/externalData/DNase/merged.DNase.pval.signal.peaks -u | wc -l >> eRNA.overlap.txt
# GWAS SNPs
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | cut -f1-3 | sort -u | intersectBed -a $tmp -b stdin -u | wc -l >> eRNA.overlap.txt

# R
setwd("~/eRNAseq")
features=c("TFBS","P300","CAGE", "Histone","VISTA", "HCNE", "DNase", "GWAS")
nHiTNE=read.table("eRNA.overlap.txt")[,1]  # c(4148,3914,1266,21823,63,15987,294,317)
for(i in 1:8){
    df=read.table(paste0("randomoverlap.",features[i],".txt"))$V1
    #cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/75490,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i])$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t", file = "randomoverlap.enrichment.xls", append =T)
    cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/71469,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i], alternative =)$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t")
}

# TFBS	4778 (6.7%)	3622 (3619,3626)	0	
# P300	3674 (5.1%)	2503 (2500,2506)	0	
# CAGE	1212 (1.7%)	942 (940,944)	0	
# Histone	20505 (28.7%)	8849 (8843,8854)	0	
# VISTA	59 (0.1%)	54 (53,54)	3.007842e-91	
# HCNE	277 (0.4%)	248 (247,249)	0	
# DNase	10894 (15.2%)	7707 (7702,7712)	0	
# GWAS	4446 (6.2%)	3293 (3289,3296)	0

################################################################################################
# for in vitro/vivo test
################################################################################################
mkdir ~/eRNAseq/invitro; cd ~/eRNAseq/invitro
# random (negative) : not overlap with any exons, promoter, eRNA, and other possible enhancer regions
cat ../*/eRNA.bed ../externalData/toExclude.bed eRNA.bed ../externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed ../externalData/CAGE/permissive_enhancers.bed ../externalData/Segment/15_coreMarks_segments.E6E7E12.bed ../externalData/VISTA/hg19.tested_regions.bed ../externalData/DNase/regions_enh_merged.brain.bed ../externalData/Conservation/HCNE_hg19_danRer7_70pc_50col.bed $GENOME/Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.bed ../externalData/Segment/wgEncodeAwgSegmentationCombinedHelas3.Enh.bed ../externalData/CAGE/hg19_permissive_enhancer.HeLaS3.bed | cut -f1-3 | sortBed | mergeBed -i - | bedtools shuffle -excl stdin -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt | awk '{OFS="\t"; print $0, $3-$2}' > negative.bed

awk '{OFS="\t"; $4=$1"_"$2"_"$3"_"($3-$2); if($5>500 && $5<800) print;}' negative.bed | head -n20 | bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed - -fo stdout | sed 's/_/\t/g' > negative.xls

R
#for(i in c("HCILB_SNDA","HC_MCPY","HC_TCPY","HC_PBMC","HC_FB")){
for(i in c("HCILB_SNDA","HC_PY","HC_nonNeuron")){
  print(i)
  # add privacy and mean +/- std for each group 
  df1=read.table(paste0('~/eRNAseq/',i,'/eRNA.meanRPM.xls'), header=T)
  rownames(df1) = df1[,1]; df1=df1[,-1]
  df1=data.frame(median=round(apply(df1, 1, median),2), mean=round(apply(df1, 1, mean),2), sd=round(apply(df1, 1, sd),2))
  df2=read.table(paste0('~/eRNAseq/',i,'/eRNA.characterize.xls'), header=T)
  df3=read.table(paste0('~/eRNAseq/',i,'/eRNA.privacy.bed'), stringsAsFactors = F); colnames(df3) = c('chr','start','end','ID','private')
  dim(df1); dim(df2); dim(df3)
  # merge into a big table
  df=cbind(df3, df1[df3$ID, ], df2[df3$ID, ])
  write.table(df, paste0('~/eRNAseq/',i,'/eRNA.characterize2.xls'), sep="\t", quote = F, col.names = NA, row.names = T)
  
  # top 10 highly expressed HTNE in each group (mean +/- std)
  df=read.table(paste0('~/eRNAseq/',i,"/eRNA.meanRPM.xls"), header=T)
  rownames(df) = df[,1]; df=df[,-1]
  df=data.frame(type=i, median=round(apply(df, 1, median),2), mean=round(apply(df, 1, mean),2), sd=round(apply(df, 1, sd),2))
  df=head(df[with(df, order(-median, sd)),],10)
  write.table(df, paste0("top10highlyexpressed.HTNE.",i,'.xls'), quote =F, sep="\t")
}

################################################################################################
# eRNA cluster
################################################################################################
# distribution of distance between neighboring eRNA (in example of SNDA)
cat HCILB_SNDA/eRNA.bed | sortBed | awk '{if(chr==$1) print $2-end; end=$3;chr=$1;}' > HCILB_SNDA/eRNA.gap.length.txt
Rscript ~/pipeline/src/_HTNE.distance.threshold.R HCILB_SNDA/eRNA.gap.length.txt

# clusters with different distance gap threshold
for i in HCILB_SNDA HC_PY HC_nonNeuron HC_PBMC HC_FB HC_MCPY HC_TCPY; do
    echo $i;
    bedtools merge -d 100000 -i $i/eRNA.bed -c 4 -o count | awk '{OFS="\t"; print $0, $3-$2}' | sort -k5,5nr | intersectBed -a - -b <(awk '{OFS="\t"; if($8=="protein_coding") print $1,$2,$3,$7,$5,$6}' $GENOME/Annotation/Genes/genes.bed) -wao | groupBy -g 1-5 -c 9 -o distinct > $i/eRNA.cluster.bed
    bedtools merge -d 50000 -i $i/eRNA.bed -c 4 -o count | awk '{OFS="\t"; print $0, $3-$2}' | sort -k5,5nr | intersectBed -a - -b <(awk '{OFS="\t"; if($8=="protein_coding") print $1,$2,$3,$7,$5,$6}' $GENOME/Annotation/Genes/genes.bed) -wao | groupBy -g 1-5 -c 9 -o distinct > $i/eRNA.cluster.d50k.bed
    bedtools merge -d 30000 -i $i/eRNA.bed -c 4 -o count | awk '{OFS="\t"; print $0, $3-$2}' | sort -k5,5nr | intersectBed -a - -b <(awk '{OFS="\t"; if($8=="protein_coding") print $1,$2,$3,$7,$5,$6}' $GENOME/Annotation/Genes/genes.bed) -wao | groupBy -g 1-5 -c 9 -o distinct > $i/eRNA.cluster.d30k.bed
done

# venn diagram for genes in the top 100 clusters 
for i in HCILB_SNDA HC_PY HC_nonNeuron; do head -n100 $i/eRNA.cluster.d50k.bed | cut -f6 | sed 's/,/\n/g' > $i/eRNA.cluster.d50k.bed.top100.genelist; done
R
library(VennDiagram) # install.packages('VennDiagram')
setwd("~/eRNAseq")
a=read.table("HCILB_SNDA/eRNA.cluster.d50k.bed", header=F, stringsAsFactors = F, nrows =100)
b=read.table("HC_PY/eRNA.cluster.d50k.bed", header=F, stringsAsFactors = F, nrows =100)
c=read.table("HC_nonNeuron/eRNA.cluster.d50k.bed", header=F, stringsAsFactors = F, nrows =100)

overlaps = calculate.overlap(list(SNDA = unlist(strsplit(a$V6, ",")), PY = unlist(strsplit(b$V6, ",")), nonNeuron = unlist(strsplit(c$V6, ","))))
sapply(overlaps, length)
names(overlaps) = c("abc",'ab','ac','bc','a','b','c'); # <-- "a5" "a2" "a4" "a6" "a1" "a3" "a7"
## run eulerAPE_3 in command line or in client
#java -Duser.language=en -Duser.region=en -jar ~/tools/eulerAPE_3.0.0.jar -i "areas_example.els" -o "./" -l yes -c no --curves circles 

# enriched gene set for the top genes
cat(overlaps$a)  # --> http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp (C2:Canonical pathways, 1330 gene sets, top 10 with FDR<0.05) --> https://docs.google.com/spreadsheets/d/12_NJxoux3Zey-ZbSBLM-w0cdioLPPjkEw73pldh3_ek/edit#gid=0 
require(RCurl)
gsea=read.delim(textConnection(getURL("https://docs.google.com/spreadsheets/d/12_NJxoux3Zey-ZbSBLM-w0cdioLPPjkEw73pldh3_ek/pub?output=tsv")), stringsAsFactors =F)
gsea=gsea[nrow(gsea):1, ]
head(gsea)
gsea$p.value=as.numeric(gsub(" ","",gsea$p.value)); gsea$FDRq.value=as.numeric(gsub(" ","",gsea$FDRq.value));
pdf("~/eRNAseq/GSEA.HTNEcluster.host.genes.pdf", width=9, height=nrow(gsea)/4)
gsea$Description <- factor(gsea$Description, unique(as.character(gsea$Description)))
library(ggplot2)
cl=read.delim(textConnection(getURL("https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/pub?gid=1995457670&single=true&output=tsv")), stringsAsFactors =F)
head(cl)
colours=paste0("#",cl$HEX[match(c('HCILB_SNDA','HC_PY', 'HC_nonNeuron', 'HC_Neuron'), cl$ITEM)])
names(colours) = c('SNDA only','PY only', 'Non-neuron only', 'SNDA & PY')

p=ggplot(gsea, aes(x = Description, y = -log10(p.value), fill=Geneset)) + scale_fill_manual(values=colours) + geom_bar(stat = "identity") + coord_flip() + theme_bw() +  theme_classic() + theme(legend.position="none")
print(p)
dev.off() 

# run topGO for the gene
source("~/pipeline/bin/lib.R")
topGOenrichment(overlaps$a, topN=100, pCutoff=0.001, type='all', output='GOenrichment.a')
topGOenrichment(overlaps$b, topN=100, pCutoff=0.001, type='all', output='GOenrichment.b')
topGOenrichment(overlaps$c, topN=100, pCutoff=0.001, type='all', output='GOenrichment.c')
topGOenrichment(overlaps$ab, topN=100, pCutoff=0.001, type='all', output='GOenrichment.ab')

# R end

## length comparison of top100 clusters
for i in HCILB_SNDA HC_PY HC_nonNeuron; do head -n100 $i/eRNA.cluster.d50k.bed | awk -vi=$i '{OFS="\t"; print $1,$2,$3,i":"$5":"$6}'; done | sortBed | mergeBed -c 4 -o distinct -delim "|" | awk '{OFS="\t"; split($4,a,"|"); for(i=1;i<=length(a);i++) print $1,$2,$3,a[i];}' | sed 's/:/\t/g' | sort -k1,1 -k2,2 -k3,3nr | awk '{if(id!=$1$2) print; id=$1$2; }' 

for i in HCILB_SNDA HC_PY HC_nonNeuron; do head -n100 $i/eRNA.cluster.d50k.bed | awk -vi=$i '{OFS="\t"; print $1,$2,$3,$6}'; done | sortBed | mergeBed -c 4 -o distinct -delim "|" | intersectBed -a - -b HCILB_SNDA/eRNA.cluster.d50k.bed -wao | groupBy -g 1-4 -c 11 -o max | intersectBed -a - -b HC_PY/eRNA.cluster.d50k.bed -wao | groupBy -g 1-5 -c 12 -o max | intersectBed -a - -b HC_nonNeuron/eRNA.cluster.d50k.bed -wao | groupBy -g 1-6 -c 13 -o max > HTNE.top100clustersUnion.txt
for i in HCILB_SNDA HC_PY HC_nonNeuron; do head -n20 $i/eRNA.cluster.d50k.bed | awk -vi=$i '{OFS="\t"; print $1,$2,$3,$6}'; done | sortBed | mergeBed -c 4 -o distinct -delim "|" | intersectBed -a - -b HCILB_SNDA/eRNA.cluster.d50k.bed -wao | groupBy -g 1-4 -c 11 -o max | intersectBed -a - -b HC_PY/eRNA.cluster.d50k.bed -wao | groupBy -g 1-5 -c 12 -o max | intersectBed -a - -b HC_nonNeuron/eRNA.cluster.d50k.bed -wao | groupBy -g 1-6 -c 13 -o max > HTNE.top20clustersUnion.txt

R
setwd("~/eRNAseq")
df=read.table("HTNE.top20clustersUnion.txt", header=F)
#colnames(df) =c('chr','start','end','genes','SNDA','PY','Non-euron')
pdf("HTNE.top20clustersUnion.pdf", width=5, height=5)
par(mar=c(4,4,1,1))
with(df, {plot(V5,V6, cex=.8,xlab="SNDA",ylab="PY", xlim=range(V5,V6), ylim=range(V5,V6), asp=1, pch=21, col='white', bg=rgb((V5-min(V5))/(max(V5)-min(V5)),(V6-min(V6))/(max(V6)-min(V6)),(V7-min(V7))/(max(V7)-min(V7))))})
with(df, {plot(V5,V7, cex=.8,xlab="SNDA",ylab="Non-neuron", xlim=range(V5,V6), ylim=range(V5,V6),asp=1, pch=21, col='white', bg=rgb((V5-min(V5))/(max(V5)-min(V5)),(V6-min(V6))/(max(V6)-min(V6)),(V7-min(V7))/(max(V7)-min(V7))))})
with(df, {plot(V6,V7, cex=.8,xlab="PY",ylab="Non-neuron", xlim=range(V5,V6), ylim=range(V5,V6),asp=1, pch=21, col='white', bg=rgb((V5-min(V5))/(max(V5)-min(V5)),(V6-min(V6))/(max(V6)-min(V6)),(V7-min(V7))/(max(V7)-min(V7))))})
library(scatterplot3d) # install.packages('scatterplot3d')
with(df, {
  r=(V5-min(V5))/(max(V5)-min(V5));
  g=(V6-min(V6))/(max(V6)-min(V6));
  b=(V7-min(V7))/(max(V7)-min(V7))
  scatterplot3d(V5/1000,V6/1000,V7/1000, angle=50, cex.symbols=.7, xlim=range(V5/1000,V6/1000,V7/1000), ylim=range(V5/1000,V6/1000,V7/1000),  xlab="SNDA",ylab="PY",zlab="Non-neuron",pch=21, color='white', lwd=0.7, bg=rgb(r,g,b))})
with(df, {
  r=1-((V6-min(V6))/(max(V6)-min(V6)) + (V7-min(V7))/(max(V7)-min(V7)))/2;
  g=1-((V5-min(V5))/(max(V5)-min(V5)) + (V7-min(V7))/(max(V7)-min(V7)))/2;
  b=1-((V5-min(V5))/(max(V5)-min(V5)) + (V6-min(V6))/(max(V6)-min(V6)))/2;
  scatterplot3d(V5/1000,V6/1000,V7/1000, angle=50, cex.symbols=.7, xlim=range(V5/1000,V6/1000,V7/1000), ylim=range(V5/1000,V6/1000,V7/1000),  xlab="SNDA",ylab="PY",zlab="Non-neuron",pch=21, color='black', lwd=0.7, bg=rgb(r,g,b))})
dev.off()


# intronic fpkm
for i in HCILB_SNDA HC_PY HC_nonNeuron HC_PBMC HC_FB HC_MCPY HC_TCPY; do
    echo $i;
    awk '{OFS="\t"; $4=$4"="$1"_"$2"_"$3; print}' $ANNOTATION/introns.meta.bed | bigWigAverageOverBed ~/neurogen/rnaseq_PD/results/merged/trimmedmean_backup/trimmedmean.uniq.normalized.$i.bw stdin stdout 2> /dev/null | sed 's/=[^\t]*//g' | groupBy -g 1 -c 2,4 -o sum,sum | sed 's/___/\t/g' | cut -f2,4- | awk 'BEGIN{OFS="\t";print "id","intronLength","totalRPM","meanRPM"}{print $0,$3/$2}' > $i/meanRPM.of.metaintron.by.gene.tab &
done

################################################################################################
# HITNE intron co-transcriptional splicing?
################################################################################################
for i in HC_PBMC HC_FB HC_MCPY HC_TCPY HCILB_SNDA; do
    echo $i;
    intersectBed -a $ANNOTATION/introns.meta.bed -b eRNA.bed -u | awk '($3-$2)>50000' | $pipeline_path/bin/toBinRegionsOnBigwig.sh /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw - 100 > $i.HiTNE.introns.100bin.rpm &
    awk '($3-$2)>50000' $ANNOTATION/introns.meta.bed | $pipeline_path/bin/toBinRegionsOnBigwig.sh /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw - 100 > $i.all.introns.100bin.rpm &
done

R
pdf("cotranscriptionalsplicing.intron.pdf", width=6, height=4, paper='usr')
celltypes=c("HC_FB", "HC_PBMC", "HC_MCPY", "HC_TCPY","HCILB_SNDA")
colors=c("#2C2C2C", "#5F5F5F", "#99D8C9", "#9EBCDA","#2CA25F")
for(i in 1:5){
    df=read.table(paste0(celltypes[i],".HiTNE.introns.100bin.rpm"), header=F)[,-1]
    if(i==1) plot(apply(df[,1:99], 2, mean, trim=.1), frame.plot=F, las=1, type='l',lwd=3, col=colors[i], ylim=c(0,0.07), ylab="Average RPM", xlab="Meta-intron position (% of total length)") else points(apply(df[,1:99], 2, mean, trim=.1), type='l',lwd=3, lty=1, col=colors[i])
    df=read.table(paste0(celltypes[i],".all.introns.100bin.rpm"), header=F)[,-1]
    points(apply(df[,1:99], 2, mean, trim=.1), type='l',lwd=3, lty=3, col=colors[i])
}
abline(v=c(0,100),lty=2)
legend("topright",c(celltypes), col=colors, pch='-',lty=1, cex=.8, lwd=4, bty='n')
legend("top",c("HiTNE introns  ","all introns"), pch='-',lty=c(1,3), cex=.8, lwd=4, bty='n')

celltypes=c("HCILB_SNDA")
colors=c("#2CA25F")
for(i in 1){
    df=read.table(paste0(celltypes[i],".HiTNE.introns.100bin.rpm"), header=F)[,-1]
    if(i==1) plot(apply(df[,1:99], 2, mean, trim=.1), frame.plot=F, las=1, type='l',lwd=3, col=colors[i], ylim=c(0,0.02), ylab="Average RPM", xlab="Meta-intron position (% of total length)") else points(apply(df[,1:99], 2, mean, trim=.1), type='l',lwd=3, lty=1, col=colors[i])
    df=read.table(paste0(celltypes[i],".all.introns.100bin.rpm"), header=F)[,-1]
    points(apply(df[,1:99], 2, mean, trim=.1), type='l',lwd=3, lty=3, col=colors[i])
}
abline(v=c(0,100),lty=2)
legend("topright",c("HiTNE introns  ","all introns"), pch='-',lty=c(1,3), cex=.8, lwd=4, bty='n')
dev.off()


rm introns.meta.*; cat $ANNOTATION/introns.meta.bed | awk '{OFS="\t"; if($4==id){ sup=($6=="+")?(end-500):$2; sdown=($6=="+")?$2:(end-500); print $1,sup,sup+500,$1"_"$2"_"$3 >> "introns.meta.up.bed"; print $1,sdown,sdown+500,$1"_"$2"_"$3 >> "introns.meta.down.bed"; end=$3;} else {if(($3-$2)>500) {id=$4;end=$3;}}}'

for i in HC_PBMC HC_FB HC_MCPY HC_TCPY HCILB_SNDA; do
    echo $i;
    #bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw introns.meta.up.bed stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.all.introns.up500.rpm &
    #bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw introns.meta.down.bed stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.all.introns.down500.rpm &
    #cut -f2 Hostgene.length.nHITNE.metaExon.metaIntron.xls | fgrep -f - introns.meta.up.bed | bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw stdin stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.HiTNE.introns.up500.rpm &
    #cut -f2 Hostgene.length.nHITNE.metaExon.metaIntron.xls | fgrep -f - introns.meta.down.bed | bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw stdin stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.HiTNE.introns.down500.rpm &
    awk '$5>10' Hostgene.length.nHITNE.metaExon.metaIntron.xls | cut -f2 | fgrep -f - introns.meta.up.bed | bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw stdin stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.HiTNE.gt10.introns.up500.rpm &
    awk '$5>10' Hostgene.length.nHITNE.metaExon.metaIntron.xls | cut -f2 | fgrep -f - introns.meta.down.bed | bigWigAverageOverBed /data/neurogen/rnaseq_PD/results/merged/trimmedmean.uniq.normalized.$i.bw stdin stdout | cut -f1,5 | awk '{OFS="\t"; split($1,a,"_"); print $1,a[3]-a[2],$2;}'> $i.HiTNE.gt10.introns.down500.rpm &
done

# R
df=c()
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    df=rbind(df, data.frame(celltype=i,region="up500",read.table(paste0(i,".all.introns.up500.rpm"), header=F)[,c(2,3)]))
    df=rbind(df, data.frame(celltype=i,region="down500",read.table(paste0(i,".all.introns.down500.rpm"), header=F)[,c(2,3)]))
}
colnames(df)=c('celltype','region',"intronLength","RPM")
df=cbind(df, Length=cut(df$intronLength, breaks=c(0,10000,50000,100000,2000000), labels=c("S","M","L","XL")))

pdf("cotranscriptionalsplicing2.all.pdf", width=4, height=4, paper="us")
meanRPM=tapply(df$RPM, list(df$Length,df$region, df$celltype), mean)
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    barplot(t(meanRPM[,,i]),beside=T, xlab="Length of downstream intron", ylab="Average RPM", main=i, col=c(0,1), las=1)
}
dev.off()
        
df=c()
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    df=rbind(df, data.frame(celltype=i,region="up500",read.table(paste0(i,".HiTNE.introns.up500.rpm"), header=F)[,c(2,3)]))
    df=rbind(df, data.frame(celltype=i,region="down500",read.table(paste0(i,".HiTNE.introns.down500.rpm"), header=F)[,c(2,3)]))
}
colnames(df)=c('celltype','region',"intronLength","RPM")
df=cbind(df, Length=cut(df$intronLength, breaks=c(0,10000,50000,100000,2000000), labels=c("S","M","L","XL")))

pdf("cotranscriptionalsplicing2.HiTNE.pdf", width=4, height=4, paper="us")
meanRPM=tapply(df$RPM, list(df$Length,df$region, df$celltype), mean)
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    barplot(t(meanRPM[,,i]),beside=T, xlab="Length of downstream intron", ylab="Average RPM", main=i, col=c(0,1), las=1)
}
dev.off()

df=c()
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    df=rbind(df, data.frame(celltype=i,region="up500",read.table(paste0(i,".HiTNE.gt10.introns.up500.rpm"), header=F)[,c(2,3)]))
    df=rbind(df, data.frame(celltype=i,region="down500",read.table(paste0(i,".HiTNE.gt10.introns.down500.rpm"), header=F)[,c(2,3)]))
}
colnames(df)=c('celltype','region',"intronLength","RPM")
df=cbind(df, Length=cut(df$intronLength, breaks=c(0,10000,50000,100000,2000000), labels=c("S","M","L","XL")))

pdf("cotranscriptionalsplicing2.HiTNE.gt10.pdf", width=4, height=4, paper="us")
meanRPM=tapply(df$RPM, list(df$Length,df$region, df$celltype), mean)
for(i in c("HC_PBMC", "HC_FB", "HC_MCPY", "HC_TCPY", "HCILB_SNDA")){
    barplot(t(meanRPM[,,i]),beside=T, xlab="Length of downstream intron", ylab="Average RPM", main=i, col=c(0,1), las=1)
}
dev.off()        


#############################################################
# clustering of eRNA expression
#############################################################
Rscript eRNA.clustering.R eRNA.80samples.RPKM.xls
Rscript eRNA.clustering.R eRNA.80samples.rawcount.xls


#############################################################
# DESeq of eRNA 
#############################################################



#############################################################
# target genes of eRNA assigned by correlation of expression
#############################################################

cut -f1,9- ~/neurogen/rnaseq_PD/results/merged/genes.fpkm.allSamples.uniq.xls > genes.fpkm.allSamples.uniq.xls
bsub -q big-multi -n 4 Rscript /PHShome/xd010/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R eRNA.meanRPM.xls /PHShome/xd010/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls eRNA.correlate.gene.in.fpkm.rho.tab
bsub -q big-multi -n 4 Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R ~/eRNAseq/HCILB_SNDA/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/genes.fpkm.cuffnorm.allSamples.uniq.xls ~/eRNAseq/HCILB_SNDA/eRNA.correlate.gene.allSamples.rho.tab

# join pcc and rho
join -j 1 <(awk '{print $1"__"$2,$3;}' eRNA.correlate.gene.in.RPKM.pcc.tab|sort) <(awk '{print $1"__"$2,$3;}' eRNA.correlate.gene.in.RPKM.rho.tab|sort) | sort -r | sed 's/__/\t/g;s/ /\t/g' >eRNA.correlate.gene.in.RPKM.cor.tab

# eRNA per gene
cut -f2 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | datamash mean 1
#9.451491660794
# gene per eRNA
cut -f1 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sed 's/^\s*//g;s/ /\t/g' | datamash mean 1
#15.613116026387

# top eRNA with most target genes
cut -f1 eRNA.correlate.gene.in.RPKM.tab | sort | uniq -c | sort -k1,1nr

# xyplot

grep chr eRNA.correlate.gene.in.RPKM.rho.tab | head -n20 | while read x y rest ; do
echo $x, $y;
set -v
    grep -P "$x|locus" eRNA.allsamples.RPKM.tab > /tmp/xyplot.tab
    grep $y genes.fpkm.allSamples.uniq.xls >> /tmp/xyplot.tab
    Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.xyplot.R /tmp/xyplot.tab "$x.$y"
set +v
done


# RPM with p=0.001 in all samples is: 0.64
# RPM with p=0.001 in HC_SNDA.trimmedmean.uniq.normalized is: 0.17

#3: robust set of eRNA


# length distribution
for i in *.trimmedmean.uniq.normalized.eRNA.bed; do echo $i; wc -l $i; awk '{print $3-$2}' $i | textHistogram -binSize=20 -maxBinCount=50 stdin; done

awk '{OFS="\t"; if(($3-$2)>=200) print $1, $2, $3, $1"_"$2"_"$3, $4}' HC_SNDA.trimmedmean.uniq.normalized.eRNA.bed > eRNA.bed

# Total counts of CAGE reads
toBinRegionsOnBigwig.sh ../CAGE/ctssTotalCounts.fwd.bw eRNA.bed 1 max > eRNA.CAGE.fwd.bed &
toBinRegionsOnBigwig.sh ../CAGE/ctssTotalCounts.rev.bw eRNA.bed 1 max > eRNA.CAGE.rev.bed &

# TF count
toBinRegionsOnBigwig.sh ../TFBS/TFBS.bigwig eRNA.bed 1 max > eRNA.TFBS.bed &

# Histone
toBinRegionsOnBigwig.sh ../Histone/Histone.SN.H3K27ac.bigwig eRNA.bed 1 max > eRNA.SN.H3K27ac.bed &

# DNase
toBinRegionsOnBigwig.sh ../DNase/DNase.bigwig eRNA.bed 1 max > eRNA.DNase.bed &

echo -e "position\tRNAseq\tCAGE.fwd\tCAGE.rev\tDNase\tH3K27ac\tTFBS" > eRNA_merged.txt
paste eRNA.*bed | sed 's/ /\t/g' | cut -f4,5,7,9,11,13,15 >> eRNA_merged.txt

R
df=read.table("eRNA_merged.txt", header=T)
rownames(df)=df[,1]; df=df[,-1]
attach(df)
x=log10(1+cbind(CAGE.fwd, CAGE.rev))
plot(x, pch=1, col='#ff000022', cex=sqrt(RNAseq))

library(flashClust)
d=df[,grep('CAGE',colnames(df))]
d=d[sample(nrow(d),2000),]
dis=dist(d)
hc <- hclust(dis)
plot(hc)

df=df[with(df, order(RNAseq)),]
for(i in 1:ncol(df)){
    image(t(df[,i,drop=F]))
}
