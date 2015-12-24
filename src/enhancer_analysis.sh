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

# for minor cell tyes
intersectBed -a HCILB_SNDA/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HC_MCPY/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HCILB_SNDA/eRNA.private.minor.bed
intersectBed -a HC_TCPY/eRNA.bed -b <(cat HCILB_SNDA/eRNA.bed HC_MCPY/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_TCPY/eRNA.private.minor.bed
intersectBed -a HC_MCPY/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_FB/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_MCPY/eRNA.private.minor.bed
intersectBed -a HC_FB/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_MCPY/eRNA.bed HC_PBMC/eRNA.bed) -v > HC_FB/eRNA.private.minor.bed
intersectBed -a HC_PBMC/eRNA.bed -b <(cat HC_TCPY/eRNA.bed HCILB_SNDA/eRNA.bed HC_FB/eRNA.bed HC_MCPY/eRNA.bed) -v > HC_PBMC/eRNA.private.minor.bed

# shared btw SNDA and the other
> HCILB_SNDA/eRNA.shared.minor.txt
for i in HC_TCPY HC_MCPY HC_PBMC HC_FB; do echo $i; echo "HCILB_SNDA vs. $i:" `intersectBed -a HCILB_SNDA/eRNA.bed -b $i/eRNA.bed -u | wc -l` >> HCILB_SNDA/eRNA.shared.minor.txt; done


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


################################################################################################
# eQTL of eRNA
################################################################################################
cd HCILB_SNDA;
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.R # for HCILB_SNDA only
## run permutations in bash
for i in `seq 1 10000`; do echo $i; [ -e permutations/permutation$i.txt ] || bsub -n 1 -M 500 -q short -J $i Rscript ~/neurogen/pipeline/RNAseq/modules/_eQTL_permutation_minP.R $i data.RData expression.postSVA.xls 1000; done
# post-eQTL analysis
bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.eQTL.after.R

## eQTL filters:
# 1. GWAS disease associated
snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
[ -e $snps_in_LD.autosomal.associations.bed ] || awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sortBed | uniq > $snps_in_LD.autosomal.associations.bed
more +2 final.cis.eQTL.xls | awk '{OFS="\t"; split($2,a,"_"); print a[1],$7-1,$7,$1,$2,$3,$4,$5;}' | intersectBed -a - -b $snps_in_LD.autosomal.associations.bed -wo | cut -f1-8,12 > final.cis.eQTL.GWAS.xls
## out of 495085 GWAS associated SNPs, 440565 (89%) are covered by our SNP sets (genotyping + imputation)

# 2. significant after multi-test correction (FDR or Bonferroni)
R
df=read.table("final.cis.eQTL.GWAS.xls", header=F, stringsAsFactors =F, sep="\t", quote="\"")
colnames(df) = c('chr','start','end','SNP','gene','beta','t.stat','p.value','trait')
df$FDR=p.adjust(df$p.value,method='fdr'); df$bonferroni=p.adjust(df$p.value,method='bonferroni')
write.table(df, "final.cis.eQTL.GWAS.adjusted.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)
quit('no')
# 3. disrupt the TFBS motif
TFBS=../externalData/TFBS/factorbookMotifPos.v3.bed
more +2 final.cis.eQTL.GWAS.adjusted.xls | intersectBed -a - -b $TFBS -wo | awk -vFS="\t" '$11<=0.05' > final.cis.eQTL.GWAS.Bonferroni5pt.TFBS.xls
more +2 final.cis.eQTL.GWAS.adjusted.xls | intersectBed -a - -b $TFBS -wo | awk -vFS="\t" '$10<=0.05' > final.cis.eQTL.GWAS.FDR5pt.TFBS.xls

cat final.cis.eQTL.GWAS.Bonferroni5pt.TFBS.xls | cut -f9 | sort | uniq -c | sort -k1,1n
grep Parkinson final.cis.eQTL.GWAS.Bonferroni5pt.TFBS.xls




awk '{OFS="\t"; split($2,a,"_"); if($7>=a[2] && $7<=a[3]) print;}' final.cis.eQTL.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<2 || $9>10) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo
more +2 final.cis.eQTL.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05 && $7>=a[2] && $7<=a[3]) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo
more +2 final.cis.eQTL.xls | sort -k2,2 | join -a 1 -1 2 -2 1 -o "1.1,0,1.3,1.4,1.5,1.6,1.7,2.28,2.6" - <(sort eRNA.characterize.xls) | sed 's/ /\t/g' | awk '{OFS="\t"; if($8<3 || $9>10) {split($2,a,"_"); print a[1],$7-1,$7,$1,$2;}}' | intersectBed -a - -b ../externalData/TFBS/factorbookMotifPos.v3.bed -wo

## disease with significant cooccurance of HTNE-eQTL and GWAS 
more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -u  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.with.eSNP.txt
more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '{OFS="\t"; split($2,a,"_"); if($6<=0.05) print a[1],$7-1,$7,$1,$2;}' | intersectBed -a $snps_in_LD.autosomal.associations.bed -b - -v  | sort -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}'> dSNP.without.eSNP.txt

R
N=as.numeric(system("wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '", intern=T)) ## total SNPs in dbSNP
n=as.numeric(system("more +2 final.cis.eQTL.new.d1e6.p1e-2.xls | awk '$6<=0.05' | cut -f1,7 | sort -u | wc -l", intern=T))  ## total eSNP
df1=read.table("dSNP.with.eSNP.txt", header=F); rownames(df1)=df1[,1]
df2=read.table("dSNP.without.eSNP.txt", header=F); rownames(df2)=df2[,1]
df=cbind(df1, df2[rownames(df1),2]); df=df[,-1]; colnames(df)=c('observed','all')  ## only the disease with both eSNP and dSNP
results = cbind(Disease_or_Trait=rownames(df), 
                df, 
                pvalue=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$p.value), 
                OR=apply(df, 1, function(x) fisher.test(matrix(c(x[1],n-x[1], x[2], N-n-x[2]), nrow = 2), alternative='greater')$estimate))
write.table(results, "dSNP.eQTL.enrichment.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)
quit('no')



################################################################################################
# characterize eRNA
################################################################################################
# see eRNA.characterize.sh
for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY ILB_SNDA PD_SNDA; do bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.characterize.sh $i; done

for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC HC_MCPY HC_TCPY ILB_SNDA PD_SNDA; do cd $i; Rscript $pipeline_path/src/eRNA.characterize.merge.R `ls eRNA.f*.txt`; cd -; done


# if the overlap is significant or not?
for i in `seq 1 1000`;
do
    echo $i;
    bsub -J random_overlap -oo _random_overlap.$i.log -eo _random_overlap.$i.log -q normal -n 1 -M 500 $pipeline_path/src/overlap.with.random.sh
done

tmp=eRNA.bed
# TFBS hotspot
awk '$6>1' externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -c | awk '$5>5' | wc -l > eRNA.overlap.txt
# P300
awk '$4=="EP300"' externalData/TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed | intersectBed -a $tmp -b stdin -u | wc -l >> eRNA.overlap.txt
# FANTOMF CAGE
intersectBed -a $tmp -b externalData/CAGE/permissive_enhancers.bed -u | wc -l >> eRNA.overlap.txt
# Roadmpa Histone
intersectBed -a $tmp -b externalData/Segment/15_coreMarks_segments.E6E7E12.bed -u | wc -l >> eRNA.overlap.txt
# VISTA positive
grep -v NULL externalData/VISTA/hg19.tested_regions.bed | intersectBed -a $tmp -b stdin -u | wc -l >> eRNA.overlap.txt
# HCNE
intersectBed -a $tmp -b externalData/Conservation/HCNE_hg19_danRer7_70pc_50col.bed -u | wc -l >> eRNA.overlap.txt
# DNase peak from Roadmap
intersectBed -a $tmp -b externalData/DNase/merged.DNase.pval.signal.peaks -u | wc -l >> eRNA.overlap.txt
# GWAS SNPs
snps_in_LD=/data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.snps_in_LD.r2.0.4.bed
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

################################################################################################
# for in vitro/vivo test
################################################################################################
# random (negative) : not overlap with any exons, promoter, eRNA, and other possible enhancer regions
cat ../toExclude.bed eRNA.bed ../TFBS/wgEncodeRegTfbsClusteredWithCellsV3.bed ../CAGE/permissive_enhancers.bed ../Segment/15_coreMarks_segments.E6E7E12.bed ../VISTA/hg19.tested_regions.bed ../DNase/regions_enh_merged.brain.bed ../Conservation/HCNE_hg19_danRer7_70pc_50col.bed $GENOME/Annotation/GWASCatalog/gwascatalog2015Apr.gwas-clean-hg19.uniq.bed ../Segment/wgEncodeAwgSegmentationCombinedHelas3.Enh.bed ../CAGE/hg19_permissive_enhancer.HeLaS3.bed | cut -f1-3 | sortBed | mergeBed -i - | bedtools shuffle -excl stdin -noOverlapping -i eRNA.bed -g $ANNOTATION/ChromInfo.txt | awk '{OFS="\t"; print $0, $3-$2}' > negative.bed

awk '{OFS="\t"; $4=$1"_"$2"_"$3"_"($3-$2); if($5>500 && $5<800) print;}' negative.bed | head -n20 | bedtools getfasta -name -tab -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -bed - -fo stdout | sed 's/_/\t/g' > negative.xls

################################################################################################
# eRNA cluster
################################################################################################
# top 10 clusters with distance gap <100k bp
bedtools merge -d 100000 -i eRNA.bed -c 4 -o count | awk '{OFS="\t"; print $0, $3-$2}' | sort -k5,5nr | head -n10

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
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R eRNA.allsamples.RPKM.tab genes.fpkm.allSamples.uniq.xls eRNA.correlate.gene.in.RPKM
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.target.correlation.R eRNA.allsamples.readscount.tab ~/neurogen/rnaseq_PD/results/DE_DESeq2/PDvsHC/htseqcount.raw.allsamples.xls eRNA.correlate.gene.in.readscount.tab

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
