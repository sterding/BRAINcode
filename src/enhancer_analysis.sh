## script to analyse eRNA defined by HCILB RNA-seq

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

cd ~/projects/PD/results/eRNA

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
# reads mapped to eRNA
################################################################################################
# reads mapped to eRNA
cd ~/eRNAseq
for i in HCILB_SNDA HC_nonNeuron HC_PY; do echo $i; for j in `cat ~/neurogen/rnaseq_PD/results/merged/samplelist.$i`; do jj=~/neurogen/rnaseq_PD/run_output/$j/uniq/accepted_hits.bam.non-rRNA-mt.bam; echo $jj; bsub -q medium -n 2 -R 'rusage[mem=6000]' $HOME/neurogen/pipeline/RNAseq/modules/_get_readscount_per_region.sh $jj $i/eRNA.bed $jj.eRNA.rawcount; done; done
for i in ~/neurogen/rnaseq_PD/run_output/*/uniq/*non-rRNA-mt.bam.eRNA.rawcount; do j=${i/\/uniq*/}; echo ${j/*output\//} `cat $i | datamash sum 5`; done| sort | sed 's/ /\t/' > reads.in.TNEs.stat.txt


################################################################################################
# define shared and private eRNA
################################################################################################
# ven diagram for overlap of 3 major groups (a: SNDA; b:PY; c:NN)
# echo a `intersectBed -a HCILB_SNDA/eRNA.fdr.bed -b <(cat HC_PY/eRNA.fdr.bed HC_nonNeuron/eRNA.fdr.bed) -v | wc -l` > venn.diagram.fdr.txt
# echo b `intersectBed -a HC_PY/eRNA.fdr.bed -b <(cat HCILB_SNDA/eRNA.fdr.bed HC_nonNeuron/eRNA.fdr.bed) -v | wc -l` >> venn.diagram.fdr.txt
# echo c `intersectBed -a HC_nonNeuron/eRNA.fdr.bed -b <(cat HC_PY/eRNA.fdr.bed HCILB_SNDA/eRNA.fdr.bed) -v | wc -l` >> venn.diagram.fdr.txt
# echo ab `intersectBed -a HCILB_SNDA/eRNA.fdr.bed -b HC_PY/eRNA.fdr.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.fdr.bed -v | wc -l` >> venn.diagram.fdr.txt
# echo ac `intersectBed -a HCILB_SNDA/eRNA.fdr.bed -b HC_nonNeuron/eRNA.fdr.bed -u | intersectBed -a stdin -b HC_PY/eRNA.fdr.bed -v | wc -l` >> venn.diagram.fdr.txt
# echo bc `intersectBed -a HC_PY/eRNA.fdr.bed -b HC_nonNeuron/eRNA.fdr.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.fdr.bed -v | wc -l` >> venn.diagram.fdr.txt
# echo abc `intersectBed -a HCILB_SNDA/eRNA.fdr.bed -b HC_nonNeuron/eRNA.fdr.bed -u | intersectBed -a stdin -b HC_PY/eRNA.fdr.bed -u | wc -l` >> venn.diagram.fdr.txt

echo a `intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_PY/eRNA.bed HC_nonNeuron/eRNA.bed) -v | wc -l` > venn.diagram.bonferroni.txt
echo b `intersectBed -a HC_PY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_nonNeuron/eRNA.bed) -v | wc -l` >> venn.diagram.bonferroni.txt
echo c `intersectBed -a HC_nonNeuron/eRNA.bed -b <(cat HC_PY/eRNA.bed HCILB_SNDA/eRNA.bed) -v | wc -l` >> venn.diagram.bonferroni.txt
echo ab `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_PY/eRNA.bed -u | intersectBed -a stdin -b HC_nonNeuron/eRNA.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo ac `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo bc `intersectBed -a HC_PY/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HCILB_SNDA/eRNA.bed -v | wc -l` >> venn.diagram.bonferroni.txt
echo abc `intersectBed -a HCILB_SNDA/eRNA.bed -b HC_nonNeuron/eRNA.bed -u | intersectBed -a stdin -b HC_PY/eRNA.bed -u | wc -l` >> venn.diagram.bonferroni.txt
# run eulerAPE_3 client and make the figure there (http://www.eulerdiagrams.org/eulerAPE/)

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

## merge private HTNEs into one file
grep only ~/eRNAseq/*/eRNA.privacy.bed | cut -f2 -d":" > ~/eRNAseq/eRNA.private.major.merged.bed
cat ~/eRNAseq/*/eRNA.meanRPM.allSamples.xls | head -n1 >  ~/eRNAseq/eRNA.private.major.merged.meanRPM.allSamples.xls
cat ~/eRNAseq/*/eRNA.meanRPM.allSamples.xls | fgrep -w -f <(cut -f4 ~/eRNAseq/eRNA.private.major.merged.bed) >> ~/eRNAseq/eRNA.private.major.merged.meanRPM.allSamples.xls
# create heatmpa and PCA clustering for the private eRNAs merged from three major cell types
Rscript $pipeline_path/src/eRNA.private.heatmap.R ~/eRNAseq/eRNA.private.major.merged.bed ~/eRNAseq/eRNA.private.major.merged.meanRPM.allSamples.xls


## fold change and p-value for private HTNEs
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HCILB_SNDA/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HC_PY/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron
bsub -q short -n 1 Rscript $pipeline_path/src/eRNA.private.stat.R ~/eRNAseq/HC_nonNeuron/eRNA.meanRPM.allSamples.xls ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_nonNeuron ~/neurogen/rnaseq_PD/results/merged/samplelist.HC_PY ~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA

# combine characterize and ttest into one file:
cut -f4 eRNA.private.major.bed | fgrep -f - eRNA.meanRPM.allSamples.ttest.txt | sort -k1,1 | join -1 1 -2 1 -a 1 - <(grep SNDAonly eRNA.characterize2.xls | sort -k1,1) | sed 's/ /\t/g' > eRNA.private.major.foldchange.characterize.xls

cut -f4 eRNA.private.major.bed | fgrep -f - eRNA.meanRPM.allSamples.ttest.txt | awk '$2>=2 && $3<0.01 && $4>=2 && $5<0.01' | cut -f1 | fgrep -f - eRNA.characterize.xls | awk '$28==1 && $6>=10 && $7>0 && $8>0' 
# chr1_66277655_66278270	19766	90.8832	0.133147	0.0339962	36	1	1	1	2	NA	0.622140.185885	NA	0	0	0	PDE4B___ENSG00000184588.13___protein_coding	223	582063	573422	1	1	1
# chr18_13566745_13567739	349746	93.3733	0.170893	0.305512	50	2	2	2	10	NA	1.23605-0.0844749	NA	0	0	0	LDLRAD4___ENSG00000168675.14___protein_coding	144	435258	418499	2	1
# chr8_41570144_41570663	183877	65.0584	0.0954287	0.505082	32	2	1	0	9	NA	40.64260.390697	NA	0	0	0	ANK1___ENSG00000029534.15___protein_coding	63	243542	232743	0	0	1
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

## Do private genes tend to contain private TNEs?
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly | cut -f1-3 | intersectBed -a - -b <(grep SNDAonly eRNA.privacy.bed) -u | wc -l # 26
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly | cut -f1-3 | intersectBed -a - -b <(grep SNDAonly eRNA.privacy.bed) -v | wc -l #889
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly -v | cut -f1-3 | intersectBed -a - -b <(grep SNDAonly eRNA.privacy.bed) -v | wc -l #26834
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly -v | cut -f1-3 | intersectBed -a - -b <(grep SNDAonly eRNA.privacy.bed) | wc -l # 20833
# R
fisher.test(matrix(c(26,20833,889,26834),nrow=2,byrow=T),alternative='greater')  # p-value =1

## Do private TNEs tend to locate in private genes?
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly | cut -f1-3 | intersectBed -b - -a <(grep SNDAonly eRNA.privacy.bed) -u | wc -l # 133
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly | cut -f1-3 | intersectBed -b - -a <(grep SNDAonly eRNA.privacy.bed) -v | wc -l #40713
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly -v | cut -f1-3 | intersectBed -b - -a <(grep SNDAonly eRNA.privacy.bed) -v | wc -l #20993
awk '{FS=","; OFS="\t"; if(NR>1) print $10,$11,$12,$1,$18,$21}' ~/Dropbox/PDBrainMap/data/privacyDataBraincode.csv | sed 's/"//g' | grep SNDAonly -v | cut -f1-3 | intersectBed -b - -a <(grep SNDAonly eRNA.privacy.bed) | wc -l # 20833
# R
fisher.test(matrix(c(133,20833,40713,20993),nrow=2,byrow=T),alternative='greater')   # p-value =1



#############################################################
# GWAS SNP enrichment
#############################################################
# Using SNAP in the final main figure
#for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh PLINK $i; done 
for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 
for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; Rscript $pipeline_path/src/eRNA.SNP.enrichment.R SNAP $i; done 
Rscript $pipeline_path/src/_eRNA.SNP.enrichment.merge.R eRNA.SNP.full.SNAP HCILB_SNDA HC_PY HC_nonNeuron
Rscript $pipeline_path/src/_eRNA.SNP.enrichment.merge.R eRNA.SNP.enrichments.SNAP HCILB_SNDA HC_PY HC_nonNeuron

# class I+II , I, II, and III separately (Re Clemens's email on 08/23/2017)
Rscript $pipeline_path/src/eRNA.SNP.enrichment.R SNAP HCILB_SNDA ## generate a eRNA.SNP.enrichment.SNAP.HTNE123.pdf file


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
grep Parkinson $snps_in_LD.autosomal.associations.bed | intersectBed -b HCILB_SNDA/eRNA.bed -a - -wo | cut -f1-4,8 | intersectBed -a - -b <(awk '{OFS="\t"; if(NR>1) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2,$5,$6}}' HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls) -wao |  sed 's/ /__/g' | sort -k1,1 -k2,2 -k11,11g | awk '{if(s!=$2){print; s=$2}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wao


################################################################################################
# eQTL of eRNA
################################################################################################
# Note: eQTL for HTNEs was ran based on 83 samples (accidently, as one sample was missed to include). 
# We re-ran eQTL for HTNEs with 84 samples, but traits like ADHD is not included in the RTC result. 
# So, for the RTC table, we still use result for 83 samples, where in the boxplot we showed expression for 84 samples. 
# This discrepency only applies to HTNE eQTL, not gene eQTL.

# Note on 2017/03/03: Actually ADHD on RTC for 83 samples is due to a RTC bug (see _RTC.R code); its RTC should be 0.15, not 0.85.
# After correcting this, ADHD is not present in both 83 and 84 samples. 
# So, we just use 84 samples in the final result. No more discrepency!

cd HCILB_SNDA;
bsub -q big -n 2 -R 'rusage[mem=10000]' -eo eQTL.run.log -oo eQTL.run.log Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.R # for HCILB_SNDA only
ln -fs final.cis.eQTL.xls final.cis.eQTL.d1e6.p1.xls
cat final.cis.eQTL.xls | awk 'NR==1 || $5<=0.01' > final.cis.eQTL.d1e6.p1e-2.xls
cat final.cis.eQTL.xls | awk '$5<=0.01 && $6<=0.05' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls

# only for HTNE class I+II
awk '$5<3{print $4}' eRNA.class.bed | fgrep -f - final.cis.eQTL.d1e6.p1.xls > final.cis.eQTL.d1e6.p1.HTNE1n2.xls
#Rscript
df=read.table("~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1.HTNE1n2.xls", header=F,stringsAsFactors=F)
df$V6=p.adjust(df$V5, method="fdr")
df=df[df$V6<=0.05 & df$V5<=0.01,]
write.table(df, "~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1.HTNE1n2.p1e-2.FDRpt5.xls", quote=F, row.names=F, col.names=F, sep="\t")


## group significant eQTLs by gene and annotate with OMIM (Table S12)
Rscript $pipeline_path/modules/_eQTL_by_gene_annotation.R final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls TNE

## RTC 
bsub -q big-multi -n 4 -M 5000 -oo RTC.run.log -eo RTC.run.log Rscript $pipeline_path/modules/_RTC.R ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt expression.postSVA.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls HTNE
## annotate: add hostgene_GWAS_SNP and hostgene_eQTL_SNP
sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.xls | awk '{OFS="\t"; split($9,a,"_"); if(NR>1) print a[1],$10-1,$10,$0}'  | intersectBed -a - -b <(cat $GENOME/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-14,18 | sort | groupBy -g 1-11 -c 12 -o distinct | awk '{OFS="\t"; split($9,a,"_"); print a[1],$11-1,$11,$0}' | intersectBed -a - -b <(cat $GENOME/Annotation/Genes/genes.bed | awk '{OFS="\t"; print $1,$2,$3,$7,$5,$6}') -wao | cut -f4-15,19 | sort | groupBy -g 1-12 -c 13 -o distinct | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.annotated.xls


## filter rules (for display purpose only): 
# 1. for each gene, take the eSNPs with the best RTC score (if there are multiple eSNPs in LD) per GWAS SNP. 
# 2. If more than one best hits and one of them is GWAS, report the GWAS one, otherwise the one with the best eQTL pvalue (if multiple, pick one arbitrarily)
# 3. Finally, for each eSNP-GWAS pair, take the gene/TNE with the best eQTL pvalue (if multiple, pick one arbitrarily)
sed 's/ /___/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.annotated.xls | sort -k2,2r -k5,5 -k7,7gr | awk '{if(id!=$2$5) {print;rtc=$7;} else if($7==rtc) print; id=$2$5;}' |  awk '$7>=0.85' | awk '{if(id!=$2$5$6) {if(a!="") print a; if($1==$5) {print; p=0;a="";} else {p=$3;a=$0;}} else {if($1==$5) {print; p=0;a="";} if($3<p) {p=$3;a=$0;}} id=$2$5$6;}END{if(a!="") print a;}' | sort -k1,1 -k5,5 -k6,6 -k3,3g | awk '{if(id!=$1$5$6) print; id=$1$5$6}' | sed 's/___/ /g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls

## manhatten plot w/ RTC (Fig.6)
Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls all

## boxplot of top eQTL w/ RTC
cat final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls | awk '{FS="\t"; OFS="\t"; split($9,a,"_"); if(NR>1) print a[1],$11-1,$11,$1,$2}' | sortBed | intersectBed -a - -b ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed -wo| cut -f5,9 > RTC.gene.snp.list
Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_boxplot.R RTC.gene.snp.list

# simplfied table with LD infor (take the Associated_transcript (i.e TNE) with the smallest eQTL p-value per SNP/Associated_transcript_hostgene/Trait)
# Note: we use LD block called by "PLINK --blocks", which use hyplotypeviewer behind and different from pairwise LD caller like SNAP)
# Note2: in the final RTC Table S10, we defined LD block for MAPT locus using a wider cutoff of "plink --blocks --ld-window-kb 1000". See $GENOME/Annotation/Variation/1000G/chr17/README for detail.
echo "#SNP OmniID Ref:Alt Chr eSNP_host_gene Associated_transcript_hostgene Associated_transcript minP Trait RTC GWAS_SNP GWAS_SNP_pos GWAS_SNP_pvalue LD" | sed 's/ /\t/g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.simplified.xls
sed 's/ /_/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls | awk '{OFS="\t"; if($7>=0.85) {split($9,a,"_");print a[1],$11-1,$11,$1,a[1],$13,($8=="NA")?"(intergenic)":$8,$9,$3,$6,$7,$5,$10}}' | intersectBed -a - -b ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed -wo | cut -f4-13,17 | awk '{OFS="\t"; $1=$1"|"$11; print}' | cut -f1-10 | sort -k1,1 -k3,3 -k4,4 -k7,7 -k6,6g | awk '{if(id!=$1$3$4$7) {print; id=$1$3$4$7;}}' | while read rs chr rest; do chr=${chr/chr/}; rs0=${rs/\|*/}; ld=`fgrep -w $rs0 $GENOME/Annotation/Variation/1000G/LDblock/Chr$chr.LD.blocks.det | head -n1 | awk '{print $1"_"$2"_"$3}'`; echo $rs $chr $rest chr$ld; done | sed 's/ /\t/g' | awk '{OFS="\t"; print "chr"$2,$10-1,$10,$9,$0}' | intersectBed -a - -b <(sed 's/ /_/g' ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pvalue.bed) -wo | awk '$11==$22' | cut -f5-15,20 | awk '{OFS="\t"; split($1,a,"[|_]"); print a[1],a[2],a[3],$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$11}' >> final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.simplified.xls
# output all pairs with RTC>=0.85
echo "#SNP OmniID Ref:Alt Chr eSNP_host_gene Associated_transcript_hostgene Associated_transcripts minP Trait RTC GWAS_SNP GWAS_SNP_pos GWAS_SNP_pvalue LD" | sed 's/ /\t/g' > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.grouped.xls
sed 's/ /_/g' final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls | awk '{OFS="\t"; if($7>=0.85) {split($9,a,"_");print a[1],$11-1,$11,$1,a[1],$13,($8=="NA")?"(intergenic)":$8,$9,$3,$6,$7,$5,$10}}' | intersectBed -a - -b ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID.bed -wo | cut -f4-13,17 | awk '{OFS="\t"; $1=$1"|"$11; print}' | cut -f1-10 | sort -k1,1 -k3,3 -k4,4 -k7,7 -k6,6g | groupBy -g 1,2,3,4,7,8,9,10 -c 5,6 -o collapse,min | awk '{OFS="\t"; print $1,$2,$3,$4,$9,$10,$5,$6,$7,$8;}' | while read rs chr rest; do chr=${chr/chr/}; rs0=${rs/\|*/}; ld=`fgrep -w $rs0 $GENOME/Annotation/Variation/1000G/LDblock/Chr$chr.LD.blocks.det | head -n1 | awk '{print $1"_"$2"_"$3}'`; echo $rs $chr $rest chr$ld; done | sed 's/ /\t/g' | awk '{OFS="\t"; print "chr"$2,$10-1,$10,$9,$0}' | intersectBed -a - -b <(sed 's/ /_/g' ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pvalue.bed) -wo | awk '$11==$22' | cut -f5-15,20 | awk '{OFS="\t"; split($1,a,"[|_]"); print a[1],a[2],a[3],$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$11}' >> final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls.grouped.xls

## run permutations in bash
mkdir permutations; 
for i in `seq 1 10000`; do [ -e permutations/permutation$i.txt ] || bsub -n 1 -M 500 -q short -J $i Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_permutation_minP.R $i data.RData expression.postSVA.xls; done

## Top eQTL at 1M bp stepping
awk '{OFS="\t"; if($5<=1e-6) print $8,$7-1,$7,$1"|"$2,-log($5)/log(10),-log($6)/log(10);}' ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls | sortBed | intersectBed -a - -b $GENOME/Annotation/Genes/genes.bed -wo | cut -f1-2,4-6,13-14 | awk '{OFS="\t"; print int($2/1000000),$0}' | sort -k2,2 -k1,1n -k5,5gr | awk '{OFS="\t"; if(id!=$1$2) {id=$1$2;p=$5;print;} else if($5==p) print;}' | sed 's/|/\t/g' | sort -k1,1 -k2,2 -k5,5 -k6,6 -k8,8 | groupBy -g 1,2,5-9 -c 3,3,4 -o count,distinct,distinct -delim "|" | sort -k3,3 | join -1 3 -2 1 - <(cut -f1,28 ~/eRNAseq/HCILB_SNDA/eRNA.characterize2.xls | sort -k1,1) -a 1 | awk '{OFS="\t"; print "HTNE",$2,$3,$1,$11,$4,$5,$6,$7,$8,$9,$10}' | sort -k3,3 -k2,2n -k6,6gr > ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt
#cat ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt ~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt > ~/Dropbox/PDBrainMap/figures/eQTL/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.xls
## PRE: add host gene symbols for top eQTL
#awk '{OFS="\t"; if($5<=1e-6) print $8,$7-1,$7,$1"|"$2,-log($5)/log(10),-log($6)/log(10);}' ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls | sort -k1,1 -k5,5gr | intersectBed -a - -b $GENOME/Annotation/Genes/genes.bed -wo | cut -f1-5,13-14 > final.cis.eQTL.d1e6.p1e-6.SNPhostgene.txt
awk '{OFS="\t"; if($5<=1e-6) print $8,$7-1,$7,$1"|"$2,-log($5)/log(10),-log($6)/log(10);}' ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls |  sed 's/|/\t/g' | sort -k5,5 | join -1 5 -2 1 - <(cut -f1,28 ~/eRNAseq/HCILB_SNDA/eRNA.characterize2.xls | sort -k1,1) -a 1 | awk '{OFS="\t"; split($8,a,"___");if(a[1]=="NA") a[3]="NA";  print $2,$3,$4,$5"|"$1,$6,a[1],a[3]}' > ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.topSNPassociatedgene.xls

## manhanten plot for eQTL
Rscript $pipeline_path/modules/_eQTL_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls all

## boxplot of top eQTL
cut -f4,12 ~/eRNAseq/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-6.SNPhostgene.simplified.txt | cut -f1 -d"|" > topeQTL.gene.snp.list
Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_boxplot.R topeQTL.gene.snp.list

## for rs1684902 SNP (E-box of USF1 binding site)
echo -e "chr10_61632545_61633312\trs1684902:61633127:A:G_A:G" > ASE.chr10_61632545_61633312.rs1684902.list
Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_boxplot.R ASE.chr10_61632545_61633312.rs1684902.list  # output 


## boxplot of rs17649553:43994648:C:T_C:T vs. ENSG00000186868.11 (MAPT), ENSG00000120071.8 (KANSL1), ENSG00000214425.2 (LRRC37A4P) and chr17_44218414_44219042  in PD samples
cd ~/neurogen/rnaseq_PD/results/eQTL/PD_SNDA
echo "ENSG00000186868.11 rs17649553:43994648:C:T_C:T
ENSG00000120071.8 rs17649553:43994648:C:T_C:T
ENSG00000214425.2 rs17649553:43994648:C:T_C:T" > MAPT.gene.snp.list
Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_boxplot.R MAPT.gene.snp.list



## eQTL analysis
# HTNE w/ eQTL
cat final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls | cut -f2 | sort | uniq -c | sort -k2,2 |awk '{OFS="\t"; print $2,$1}' | join -a 1 -1 1 -2 1 -o "0,1.2,2.20,2.28" - <(sort eRNA.characterize.xls) | sed 's/\s/\t/g' | sort -k2,2nr  ## > Table.S8
fgrep -f <(cat final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls | cut -f2 | sort -u) eRNA.characterize.xls | cut -f28 | sort | uniq -c
#  20 1
#  17 2
# 114 3
# HTNE w/o eQTL
fgrep -v -f <(cat final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls | cut -f2 | sort -u) eRNA.characterize.xls | cut -f28 | sort | uniq -c
# 11815 1
# 11773 2
# 47283 3
# Fisher' test in R
fisher.test(matrix(c(37,23588,114,47283),nrow=2,byrow=T),alternative='greater')  # p-value = 0.9924

## GO analysis for the 102 genes harboring eQTL HTNEs
cut -f2 final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls | sort -u | grep -f - eRNA.characterize2.xls | cut -f28 | grep ___ |  sed 's/___/\t/g' | cut -f1 > final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.eTNE_hostgene.txt
# then go to GSEA (http://software.broadinstitute.org/gsea/msigdb/annotate.jsp, C5 collection, top50, FDR 0.05)

## eQTL filters:
## A: colocalize --> TFBS --> GWAS
## ----------------------------------
awk '{split($2,a,"_"); if(NR==1 || ($7>(a[2]-500) && $7<(a[3]+500))) print}' final.cis.eQTL.d1e6.p1e-2.xls  > final.cis.eQTL.d1e6.p1e-2.colocalized500nt.xls
awk '{OFS="\t"; split($2,a,"_"); if(NR>1) print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' final.cis.eQTL.d1e6.p1e-2.colocalized500nt.xls | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo > final.cis.eQTL.d1e6.p1e-2.colocalized500nt.TFBS.xls
awk '{OFS="\t"; split($2,a,"_"); if(NR>1) print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' final.cis.eQTL.d1e6.p1e-2.colocalized500nt.xls | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -u | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo ## 2

## A2: colocalize --> FDR --> TFBS (N=0)
## ----------------------------------
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' final.cis.eQTL.d1e6.p1e-2.colocalized500nt.xls | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -u # NONE

## B1: GWAS --> FDR --> TFBS
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
#df=read.table("final.cis.eQTL.GWAS.FDR5pt.xls", header=F, stringsAsFactors =F, sep="\t", quote="\"")
#colnames(df) = c('chr','start','end','SNP','gene','beta','t.stat','p.value','trait','FDR','bonferroni')
q('no')

more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$10<=0.05' > final.cis.eQTL.GWAS.FDR5pt.xls
#more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$11<=0.05' # only have 36 eQTL remained after Bonferroni correction, choose to use FDR

## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.xls # N=1

# 3. disrupt the TFBS motif (N=305)
TFBS=../externalData/TFBS/factorbookMotifPos.v3.bed
intersectBed -a final.cis.eQTL.GWAS.FDR5pt.xls -b $TFBS -wo > final.cis.eQTL.GWAS.FDR5pt.TFBS.xls
## 3.2 anyone in RTC of diseases


## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls # N=0

cat final.cis.eQTL.GWAS.FDR5pt.TFBS.xls | cut -f9 | sort | uniq -c | sort -k1,1n
grep Parkinson final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

## B2: GWAS --> p<0.01 --> colocalized --> TFBS
## ----------------------------------
awk '$8<0.01' final.cis.eQTL.GWAS.adjusted.xls > final.cis.eQTL.GWAS.p1e-2.xls
Rscript $pipeline_path/modules/_eQTL_RTC_manhanttenPlot.R final.cis.eQTL.d1e6.p1e-2.xls final.cis.eQTL.d1e6.p1e-2.FDRpt5.xls.RTC.filtered.annotated.xls all
R
setwd("~/eRNAseq/HCILB_SNDA")
inputfile="final.cis.eQTL.GWAS.p1e-2.xls" # N=162,325
df=read.table(inputfile, header=F, stringsAsFactors =F, sep="\t", quote="\"")
colnames(df) = c('chr','start','end','SNP','gene','beta','t.stat','p.value','trait','FDR','bonferroni')

pdf(sub("xls","pdf",inputfile), width=12, height = 6)
# to plot the number of colocalized eQTL SNPs per disease (e.g. eSNP in HTNEs)
df$genestart=as.numeric(sub(".*_(.*)_.*","\\1",df$gene)); df$geneend=as.numeric(sub(".*_.*_(.*)","\\1",df$gene))
df0=subset(df, start>=genestart & start<=geneend)
# remove parts in "(*)"
df0$trait = gsub(" \\(.*","", df0$trait)
n_eQTL = sort(table(df0$trait), decreasing =T)
n_eQTL = n_eQTL[n_eQTL>=2]
length(n_eQTL)
par(mar=c(10,5,2,2));
barplot(n_eQTL, width=.2, space=1, las=2, cex=.6,cex.axis=.8, cex.names =1, col='#F22A7B', border='#F22A7B', ylab="Number of disease eSNPs colocalized with HTNEs")
dev.off()
quit('no')

more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$10<=0.05' > final.cis.eQTL.GWAS.FDR5pt.xls
#more +2 final.cis.eQTL.GWAS.adjusted.xls | awk -vFS="\t" '$11<=0.05' # only have 36 eQTL remained after Bonferroni correction, choose to use FDR

## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.xls # N=1

# 3. disrupt the TFBS motif (N=348)
TFBS=../externalData/TFBS/factorbookMotifPos.v3.bed
intersectBed -a final.cis.eQTL.GWAS.FDR5pt.xls -b $TFBS -wo > final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

## Note: but none of them have any SNP colocalize with eRNA, even in +/- 500bp window
awk '{OFS="\t"; split($5,a,"_"); if($3<=a[3] && $3>=a[2]) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls  # N=0
awk '{OFS="\t"; split($5,a,"_"); if($3<=(a[3]+500) && $3>=(a[2]-500)) print}' final.cis.eQTL.GWAS.FDR5pt.TFBS.xls # N=0

cat final.cis.eQTL.GWAS.FDR5pt.TFBS.xls | cut -f9 | sort | uniq -c | sort -k1,1n
grep Parkinson final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

## C: FDR > colocalize --> TFBS: N=0
## ----------------------------------
awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05 && $7>=a[2] && $7<=a[3]) print a[1],$7-1,$7,$1,$2;}' final.cis.eQTL.d1e6.p1e-2.xls  | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo  # NONE

## D: p<0.01 --> colocalize --> TFBS: N=5
## ----------------------------------
awk '{OFS="\t"; split($2,a,"_"); if($7>=a[2] && $7<=a[3]) print;}' final.cis.eQTL.d1e6.p1e-2.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -u  # N=6
awk '{OFS="\t"; split($2,a,"_"); if($7>=a[2] && $7<=a[3]) print;}' final.cis.eQTL.d1e6.p1e-2.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<2 || $9>=5) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo > final.cis.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.xls

more +2 final.cis.eQTL.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05 && $7>=a[2] && $7<=a[3]) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo
more +2 final.cis.eQTL.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<3 || $9>10) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo

## disease with significant cooccurance of HTNE-eQTL and GWAS 
## ========================================================
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
# prune GWAS SNPs
more +2 final.cis.eQTL.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.txt
more +2 final.cis.eQTL.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.txt

R
setwd("~/eRNAseq/HCILB_SNDA")
N=as.numeric(system("wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '", intern=T)) ## total SNPs in dbSNP
#N=as.numeric(system("cut -f1-3 $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed.autosomal.associations.bed | sort -u | wc -l | cut -f1 -d' '", intern=T)) ## total dSNPs in gwas+proxy
n=as.numeric(system("more +2 final.cis.eQTL.d1e6.p1e-2.xls | awk '$6<=0.05' | cut -f1,7 | sort -u | wc -l", intern=T))  ## total eSNP
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


################################################################################################
## if the overlap of HTNE with various enhancer features is significant or not?
################################################################################################

## ++++++++++
## SNDA 
## ++++++++++

cd ~/eRNAseq/HCILB_SNDA

# run the overlapping for 1000 times
rm randomoverlap.*.txt
bsub -J "random_overlap[1-1000]" -oo _random_overlap.log -eo _random_overlap.log -q normal -n 1 -M 500 $pipeline_path/src/overlap.with.random.sh

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
setwd("~/eRNAseq/HCILB_SNDA")
features=c("TFBS","P300","CAGE", "Histone","VISTA", "HCNE", "DNase", "GWAS")
nHiTNE=read.table("eRNA.overlap.txt")[,1] 
for(i in 1:8){
    df=read.table(paste0("randomoverlap.",features[i],".txt"))$V1
    #cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/75490,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i])$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t", file = "randomoverlap.enrichment.xls", append =T)
    cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/71022,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i], alternative =)$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t")
}

# TFBS	4778 (6.7%)	3622 (3619,3626)	0	
# P300	3674 (5.2%)	2503 (2500,2506)	0	
# CAGE	1212 (1.7%)	942 (940,944)	0	
# Histone	20505 (28.9%)	8849 (8843,8854)	0	
# VISTA	63 (0.1%)	58 (57,58)	2.211971e-78	
# HCNE	277 (0.4%)	248 (247,249)	0	
# DNase	10894 (15.3%)	7707 (7702,7712)	0	
# GWAS	4446 (6.3%)	3293 (3289,3296)	0

## overlap with VISTA validated enhancers
a=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | intersectBed -a eRNA.bed -b stdin -wo | grep -v NULL -c",intern = T))
b=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | intersectBed -a eRNA.bed -b stdin -wo | grep NULL -c",intern = T))
c=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | grep -v NULL -c",intern = T))
d=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | grep NULL -c",intern = T))
c(a,b,c,d)
phyper(a-1, a+b, c+d-(a+b),c, lower.tail=F) # 0.003916193
pdf("VISTA.enrichment.pdf", paper='us')
barplot(c(d/(c+d),c/(c+d),b/(a+b),a/(a+b)),col=c('gray','red','gray','red'),names.arg=paste(c('VISTA-negative','VISTA-positive','TNE-negative','TNE-positive'),"\n",c(d,c,b,a)))
dev.off()

# # overlap with class I/II TNE
# a=as.integer(system("cut -f2-5,36 eRNA.characterize2.xls | awk '$5<3' | intersectBed -a stdin -b ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed -wo | grep -v NULL -c",intern = T))
# b=as.integer(system("cut -f2-5,36 eRNA.characterize2.xls | awk '$5<3' | intersectBed -a stdin -b ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed -wo | grep NULL -c",intern = T))
# c=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | grep -v NULL -c",intern = T))
# d=as.integer(system("cat ~/eRNAseq/externalData/VISTA/hg19.tested_regions.bed | grep NULL -c",intern = T))
# c(a,b,c,d)
# phyper(a-1, a+b, c+d-(a+b),c, lower.tail=F) # 0.09716477
# # pvalue for class I only is 0.06


## ++++++++++
## SKNSH
## ++++++++++

cd ~/eRNAseq/SKNSH
## Question: The overlap of SKNSH TNE with regulatory marks from ENCODE was significantly higher than expected by chance alone with P < XXX by permutation test?

# run the overlapping for 1000 times with ENCODE defined enhancers in SK-N-SH
> randomoverlap.sknsh_ENCODE.txt
bsub -J "SKNSH_random_overlap[1-1000]" -oo _random_overlap.log -eo _random_overlap.log -q normal -n 1 -M 500 $pipeline_path/src/overlap.with.random.sknsh_ENCODE.sh 
zcat ~/eRNAseq/externalData/sknsh_ENCODE/sknsh.*gz | cut -f1-3 | sortBed | mergeBed |  intersectBed -a eRNA.bed -b - -u | wc -l > eRNA.overlap.txt

# R
setwd("~/eRNAseq/SKNSH")
features=c("sknsh_ENCODE")
nHiTNE=read.table("eRNA.overlap.txt")[,1] 
for(i in 1){
    df=read.table(paste0("randomoverlap.",features[i],".txt"))$V1
    cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/71022,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i])$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t")
}

## ++++++++++
## PFC @ BrainGVEX
## ++++++++++

cd ~/eRNAseq/BrainGVEX
## overlap of BrainGVEX TNE with ATAC-seq peaks
# if peaks called in at least 1 sample
awk '$4>=1' ~/neurogen/psychENCODE/BrainGVEX/ATAC-seq/merged.Peaks.bedGraph | intersectBed -a eRNA.bed -b - -u | wc -l
# if peaks called in at least 5 sample
awk '$4>=5' ~/neurogen/psychENCODE/BrainGVEX/ATAC-seq/merged.Peaks.bedGraph | intersectBed -a eRNA.bed -b - -u | wc -l

## Question: The overlap of BrainGVEX TNE with ATAC-seq peaks from BrainGVEX was significantly higher than expected by chance alone with P < XXX by permutation test?

# run the overlapping for 1000 times with ATAC-seq peaks in BrainGVEX
> randomoverlap.BrainGVEX_ATACseq.txt
bsub -J "BrainGVEX_ATACseq[1-1000]" -oo _random_overlap.log -eo _random_overlap.log -q short -n 1 -M 500 $pipeline_path/src/overlap.with.random.BrainGVEX_ATACseq.sh 
awk '$4>=5' ~/neurogen/psychENCODE/BrainGVEX/ATAC-seq/merged.Peaks.bedGraph | intersectBed -a eRNA.bed -b - -u | wc -l > eRNA.overlap.txt

# R
setwd("~/eRNAseq/BrainGVEX")
features=c("BrainGVEX_ATACseq")
nHiTNE=read.table("eRNA.overlap.txt")[,1] 
for(i in 1){
    df=read.table(paste0("randomoverlap.",features[i],".txt"))$V1
    cat(features[i], paste0(nHiTNE[i]," (",round(100*nHiTNE[i]/71022,1),"%)"), paste0(round(mean(df))," (", round(t.test(df, mu=nHiTNE[i])$conf.int[1]),",",round(t.test(df, mu=nHiTNE[i])$conf.int[2]),")"), t.test(df, mu=nHiTNE[i])$p.value, "\n", sep = "\t")
}

################################################################################################
# TFBS enrichment of TNE
################################################################################################
for i in HCILB_SNDA HC_nonNeuron HC_PY; do 
  cd $i; 
  # JASPAR
  #bsub -q normal -n 1 bash $pipeline_path/src/eRNA.TFBSjaspar.enrichment.v0.sh eRNA.JASPARmotifscan.txt;  ## incomplete version (fig 2d)
  #bsub -q normal -n 1 bash $pipeline_path/src/eRNA.TFBSjaspar.enrichment.v0.sh eRNA.JASPARmotifscan.2017.txt;  ## incomplete version (fig 2d)
  ## or
  cat ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/JASPARmotifsID.TFsymbol.txt | while read JASPARID TF; do bsub -q normal -n 1 -J $TF bash $pipeline_path/src/eRNA.TFBSjaspar.enrichment.core.sh $JASPARID $TF eRNA.JASPARmotifscan.2018.txt; done
  cat eRNA.JASPARmotifscan.2018.txt.MA* > eRNA.JASPARmotifscan.2018.txt; rm eRNA.JASPARmotifscan.2018.txt.MA*;
  Rscript $pipeline_path/src/eRNA.TFBSencode.enrichment.R eRNA.JASPARmotifscan.2018.txt
  
  # ENCODE
  #bash $pipeline_path/src/eRNA.TFBSencode.enrichment.sh ~/eRNAseq/externalData/TFBS/wgEncodeRegTfbsClusteredV3.bed12 eRNA.wgEncodeRegTfbsClusteredV3.txt;
  #Rscript $pipeline_path/src/eRNA.TFBSencode.enrichment.R eRNA.wgEncodeRegTfbsClusteredV3.txt;
  
  cd -
done

### Using MEME-AME and random genomic regions with GC- and length- matched
# see https://docs.google.com/document/d/1-j7bE38Qy56ztce8zey1xH2o5sM1E2smEHb0BqrbZ9E/edit#
cd ~/neurogen/TF_scan/scan_SNDA_TNE
  
bsub -q normal -J ame-all -n 1 ame --verbose 1 --oc ame2 --control Matching_backgroundg_GC_compo.fa --bgformat 1 ~/eRNAseq/HCILB_SNDA/eRNA.seq.fa ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
bsub -q normal -J ame-all -n 1 ame --verbose 1 --oc ame --control Matching_backgroundg_GC_compo.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ~/eRNAseq/HCILB_SNDA/eRNA.seq.fa ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
bsub -q normal -J ame-private -n 1 ame --verbose 1 --oc ame-private --control Matching_backgroundg_GC_compo.private.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ~/eRNAseq/HCILB_SNDA/eRNA.private.major.fa ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt
bsub -q normal -J ame-classI -n 1 ame --verbose 1 --oc ame-classI --control Matching_backgroundg_GC_compo.classI.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ~/eRNAseq/HCILB_SNDA/eRNA.classI.fa ../JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt

# R script 
pdf("~/eRNAseq/HCILB_SNDA/AME.pdf", width=25, height = 4)
#df=read.table(pipe("grep Ranksum ~/neurogen/TF_scan/scan_SNDA_TNE/ame/ame.txt"), sep=" ")[,c(6,7,8,14,21)]  # Ranksum test
df=read.table(pipe("grep motif ~/neurogen/TF_scan/scan_SNDA_TNE/ame2/ame.txt"), sep=" ", quote = "")[,c(8,9,10,14,17)]  # default Fisher's test
colnames(df)=c("Id","Symbol","Motif","pvalue","qvalue")
df$qvalue=as.numeric(sub(")","",df$qvalue))
df=subset(df, qvalue<=0.01)
write.table(df, "~/eRNAseq/HCILB_SNDA/eRNA.vs.random.JASPAR.AME.xls", sep="\t", col.names = T, row.names = F)
df1=read.table("~/eRNAseq/HCILB_SNDA/eRNA.JASPARmotifscan.2018.txt.qvalue0.01.xls", header = T)
library(tidyverse); df1 = subset(df1, type=='eRNA') %>% arrange(pvalue) %>% head(50)
sum(df1$ID %in% df$Symbol)
df1[!(df1$ID %in% df$Symbol),]



df$pvalue = df$pvalue + min(df$pvalue[df$pvalue>0])
df$qvalue = df$qvalue + min(df$qvalue[df$qvalue>0])
barplot(-log10(df$pvalue), names.arg =df$Symbol, main = "eRNA.fa", las=2, cex.names =0.4, ylab="-log10(p-value)", width = .5, space = 2)

df=read.table(pipe("grep Ranksum ~/neurogen/TF_scan/scan_SNDA_TNE/ame-private/ame.txt"), sep=" ")[,c(6,7,8,14,21)]
colnames(df)=c("ID","Symbol","Motif","pvalue","qvalue")
df$pvalue = df$pvalue + min(df$pvalue[df$pvalue>0])
df$qvalue = df$qvalue + min(df$qvalue[df$qvalue>0])
barplot(-log10(df$pvalue), names.arg =df$Symbol, main = "eRNA.private.major.fa", las=2, cex.names =0.5, ylab="-log10(p-value)", width = .5, space = 2)
abline(h = -log10(max(df$pvalue[df$qvalue<=0.01])), col="red")

df=read.table(pipe("grep Ranksum ~/neurogen/TF_scan/scan_SNDA_TNE/ame-classI/ame.txt"), sep=" ")[,c(6,7,8,14,21)]
colnames(df)=c("ID","Symbol","Motif","pvalue","qvalue")
df$pvalue = df$pvalue + min(df$pvalue[df$pvalue>0])
df$qvalue = df$qvalue + min(df$qvalue[df$qvalue>0])
barplot(-log10(df$pvalue), names.arg =df$Symbol, main = "eRNA.classI.fa", las=2, cex.names =0.4, ylab="-log10(p-value)", width = .5, space = 2)
abline(h = -log10(max(df$pvalue[df$qvalue<=0.01])), col="red")

dev.off()

# output: ~/eRNAseq/HCILB_SNDA/AME.pdf
  
### Update: using permutation 
mkdir ~/neurogen/TF_scan/scan_SNDA_TNE
cd ~/neurogen/TF_scan/scan_SNDA_TNE
# make a random background with GC and length matched, by shuffling the TNE sequences
bsub -q short -J "shuffle[1-1000]" "python ~/BiasAway/BiasAway.py m -f ~/eRNAseq/HCILB_SNDA/eRNA.seq.fa | sed 's/background_seq_.* Background sequence for //g' > eRNA.seq.shuffled.fa.\$LSB_JOBINDEX"
# scan motif for hit
#cat ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/JASPARmotifsID.TFsymbol.txt | while read JASPARID TF; do bsub -q vshort -n 1 -J "$JASPARID[1-1000]" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/$JASPARID.meme eRNA.seq.shuffled.fa.\$LSB_JOBINDEX; done
bsub -q normal -n 1 -J "scanmotif[1-1000]" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt eRNA.seq.shuffled.fa.\$LSB_JOBINDEX
bsub -q normal -n 1 -J "scan" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt eRNA.seq.fa
# merge result
> JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.shuffled.fa.txt
for i in JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.shuffled.fa.*.bed; do echo $i; cut -f4 $i|sort|uniq -c >> JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.shuffled.fa.txt; done
# eRNA.seq.fa
bsub -q normal -n 1 -J "scanmotif" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt eRNA.seq.fa
cut -f4 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.fa.bed | sort | uniq -c > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.fa.txt
### compare to random genomic background that are GC- and length-matched (from http://opossum.cisreg.ca/GC_compo/)
# upload eRNA.seq.fa to http://opossum.cisreg.ca/GC_compo/ and return Matching_backgroundg_GC_compo.fa
bsub -q normal -n 1 -J "scanmotif" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Matching_backgroundg_GC_compo.fa
cut -f4 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.fa.bed | sort | uniq -c > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.fa.txt

bsub -q normal -n 1 -J "scanmotif" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Matching_backgroundg_GC_compo.classI.fa
cut -f4 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.classI.fa.bed | sort | uniq -c > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.classI.fa.txt

bsub -q normal -n 1 -J "scanmotif" scan_motif.sh ~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt Matching_backgroundg_GC_compo.private.fa
cut -f4 JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.private.fa.bed | sort | uniq -c > JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.private.fa.txt

# statistics
R
setwd("~/neurogen/TF_scan/scan_SNDA_TNE")
tfs=read.table("~/neurogen/TF_scan/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme/JASPARmotifsID.TFsymbol.txt", col.names = c('TF','symbol'), header = F, stringsAsFactors = F)
df1=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.fa.txt", col.names = c('count','TF'), header = F, stringsAsFactors = F)
df2=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.shuffled.fa.txt", col.names = c('count','TF'), header = F, stringsAsFactors = F)
df=data.frame(do.call("rbind", apply(tfs,1,function(x) {
  tf=as.character(x[1]); sym=as.character(x[2]);
  if(tf %in% df1$TF & tf %in% df2$TF) {
    tt=t.test(x=df2$count[df2$TF==tf], mu=df1$count[df1$TF==tf], var.equal = T);
    c(tf, sym, round(tt$statistic,1), tt$p.value, tt$null.value, tt$estimate, round(tt$conf.int[1],3), round(tt$conf.int[2],3))
  }
})), stringsAsFactors=F)
colnames(df)=c("tf", "sym", "t", "p", "obs", "est", "conf_l", "conf_r")
write.table(df, "JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_eRNA.seq.shuffled.fa.txt.statistics.xls", quote = F, sep="\t", row.names = F)

df3=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt_in_Matching_backgroundg_GC_compo.fa.txt", col.names = c('count','TF'), header = F, stringsAsFactors = F)
df=data.frame(do.call("rbind", apply(tfs,1,function(x) {
  tf=as.character(x[1]); sym=as.character(x[2]);
  if(tf %in% df1$TF & tf %in% df2$TF) {
    x=matrix(c(df1$count[df1$TF==tf], sum(df1$count[df1$TF!=tf]), df3$count[df3$TF==tf], sum(df3$count[df3$TF!=tf])),nrow = 2, dimnames = list(c("TF","No TF"), c("TNE","background")));
    tt=fisher.test(x, alternative='greater');
    c(tf, sym, c(df1$count[df1$TF==tf], sum(df1$count[df1$TF!=tf]), df3$count[df3$TF==tf], sum(df3$count[df3$TF!=tf])), tt$p.value, tt$estimate)
  }
})), stringsAsFactors=F)
colnames(df)=c("tf", "sym", 'AB','AnB', 'nAB', 'nAnB', "p", "OR")
library(tidyverse)
df = df %>% mutate(p=as.numeric(p),OR=as.numeric(OR)) %>%
  filter(p<0.01/nrow(df), OR>1) %>% arrange(p)
dim(df)

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

library(threejs) # install.packages('threejs')
library(htmlwidgets)
with(df, {
    r=1-((V6-min(V6))/(max(V6)-min(V6)) + (V7-min(V7))/(max(V7)-min(V7)))/2;
    g=1-((V5-min(V5))/(max(V5)-min(V5)) + (V7-min(V7))/(max(V7)-min(V7)))/2;
    b=1-((V5-min(V5))/(max(V5)-min(V5)) + (V6-min(V6))/(max(V6)-min(V6)))/2;
    saveWidget(scatterplot3js(V5/1000,V6/1000,V7/1000, axisLabels=c("SNDA", "Non-neuron","PY"),
                   xlim=c(0,2500), ylim=c(0,2500),zlim=c(0,1000),
                   color = rgb(r,g,b), grid=T, flip.y=F,renderer="canvas"), 
               file="HTNE.top20clustersUnion.scatterplot3js.html")
    }
    )


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
Rscript eRNA.clustering.R 

#############################################################
# classification of HTNE based on regulatory annotation
#############################################################
Rscript eRNA.classification.R


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

## BRAINcode GAGE data
bigWigAverageOverBed ~/neurogen/CAGE_PDBrainMap/output_dir/HCILB_merged_accepted_hits.plus.bw eRNA.bed stdout -minMax > eRNA.bed.CAGE_HCILB_merged_accepted_hits.plus.bw.tab
bigWigAverageOverBed ~/neurogen/CAGE_PDBrainMap/output_dir/HCILB_merged_accepted_hits.minus.bw eRNA.bed stdout -minMax > eRNA.bed.CAGE_HCILB_merged_accepted_hits.minus.bw.tab
# merge to see how many TNEs having CAGE reads in both plus and minus strand
paste eRNA.bed.CAGE_HCILB_merged_accepted_hits.minus.bw.tab eRNA.bed.CAGE_HCILB_merged_accepted_hits.plus.bw.tab | cut -f1-8,11-16 | awk '$4<0 && $10>0' | wc -l
43100
paste eRNA.bed.CAGE_HCILB_merged_accepted_hits.minus.bw.tab eRNA.bed.CAGE_HCILB_merged_accepted_hits.plus.bw.tab | cut -f1-8,11-16 | awk '$4<-10 && $10>10' | wc -l
4840


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
