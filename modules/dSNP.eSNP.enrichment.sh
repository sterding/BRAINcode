# Script to run enrichment analysis of disease SNPs and eQTL SNPs in a permutation based manner
# Referecences: Towfique Raj et al. Science (2014), and Dan L. Nicolae et al. PLoS Genetics (2010), Imge Hulur et al. BMC Genomics (2015)
# Author: Xianjun Dong (xdong@rics.bwh.harvard.edu)
# Version: 1.0
# Usage: $0 dSNP.eSNP.enrichment.sh final.cis.eQTL.d1e6.p1e-2.xls

eQTL_final=$1 # eQTL_final=~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.d1e6.p1e-2.xls


dir="${eQTL_final%/*}"
mkdir $dir/dSNP.eSNP.enrichment; 

if [ ! -f $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pruned.bed ]; then
  ## Ref ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/README.txt
  ## ---------------------------
  ## step 1: download GWAS SNPs
  ## ---------------------------
  cd $GENOME/Annotation/GWASCatalog/
  curl -s https://www.ebi.ac.uk/gwas/api/search/downloads/full > gwas_catalog_v1.0-downloaded_2015-11-04.tsv
  GWAS_snps=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
  awk '{FS="\t"; OFS="\t"; if(NR>1 && $13!="") print "chr"$12,$13-1,$13,"hg38_chr"$12"_"$13"_"$22,$28,".",$8}'  gwas_catalog_v1.0-downloaded_2015-11-04.tsv | sed 's/ /|/g' | liftOver -bedPlus=4 -tab stdin hg38ToHg19.over.chain.gz stdout unmapped | sort -k1,1 -k2,2 -k7,7 -k6,6g | awk '{if(id!=$4$7) {print; id=$4$7}}' | sed 's/|/ /g' > gwas_catalog_v1.0-downloaded.hg19.bed
  
  ## ---------------------------
  ## step 2: prune to get indepedant SNPs
  ## ---------------------------
  > gwas_catalog_v1.0-downloaded.hg19.pruned.bed
  for i in `seq 1 22` X; do 
    echo chr$i; 
    # extract subset from plink bed/bim files converted from 1000G vcf, using SNP ID
    grep chr$i gwas_catalog_v1.0-downloaded.hg19.pvalue.bed | cut -f4 | cut -f4 -d"_" > /tmp/gwas_chr$i
    [ `wc -l /tmp/gwas_chr$i | cut -f1 -d' '` == "0" ] && continue;
    # plink prune with extracted subset
    plink --bfile /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/1000G/Chr$i --extract /tmp/gwas_chr$i --indep-pairwise 50 5 0.3 --allow-no-sex --out /tmp/gwas_chr$i --silent
    # subset of gwas (Note: we have no way to prune SNPs by disease. So, if a SNP associates with multiple diseases, it will be kept or removed as a whole)
    sed 's/_/|/g' gwas_catalog_v1.0-downloaded.hg19.pvalue.bed | fgrep -w -f /tmp/gwas_chr$i.prune.in | sed 's/|/_/g' >> gwas_catalog_v1.0-downloaded.hg19.pruned.bed
  done
fi

## ---------------------------
## step 3: compute MAF, number of LD partners (r2>0.8), and distance from transcription start site (TSS) for all genotyped SNPs
## ---------------------------
SNPpos=/data/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID
cd $dir/dSNP.eSNP.enrichment

if [ ! -f ALL.combined ]; then
  # all SNPs on the array
  [ -f $SNPpos.bed ] || awk '{OFS="\t"; print $2,$3-1,$3,$1}' $SNPpos > $SNPpos.bed
  #for i in 1 `seq 1 22`; do echo $i; bsub -q normal -n 1 bash $pipeline_path/src/_SNP.maf.ld.d2tss.sh $SNPpos.bed $i ALL_chr; done
  bsub -J "maf[1-22]" -q normal -n 1 bash $pipeline_path/src/_SNP.maf.ld.d2tss.sh $SNPpos.bed ALL_chr
  [ ] && cat ALL_chr*.combined | awk '{OFS="\t"; $2=int($2/0.05); $3=($3>10)?10:$3;$4=($4>200000)?200000:$4; $4=int($4/10000); print $1, "MAF"$2"LD"$3"TSS"$4;}' > ALL.combined
  awk '{print $1 >> $2}' ALL.combined
fi

## ---------------------------
## step 4: calculate number of eSNPs that are also GWAS dSNPs (N_deSNP)
## ---------------------------

# use FDR<0.05 as cutoff to define eSNP
awk '$6<=0.05' $eQTL_final | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' > ${eQTL_final/xls/FDR5e-2.bed}
# # pvalue < 0.001
# awk '$5<=0.001' $eQTL_final | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' > ${eQTL_final/xls/FDR1e-3.bed}
# # pvalue < 0.005
# awk '$5<=0.005' $eQTL_final | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' > ${eQTL_final/xls/FDR5e-3.bed}
intersectBed -a ${eQTL_final/xls/FDR5e-2.bed} -b $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pruned.bed -wo | cut -f12 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > N_deSNP.txt

## ---------------------------
## step 5: generate 1000 random sets of SNPs with matched MAF, number of LD partners (r2>0.8), and distance from transcription start site (TSS) with eSNP. For each set, calculate number of each eSNP set that are also dSNPs (n_deSNP)
## ---------------------------
# only eSNPs
bsub -J "maf[1-22]" -q normal -n 1 bash $pipeline_path/src/_SNP.maf.ld.d2tss.sh ${eQTL_final/xls/FDR5e-2.bed} esnp_chr

cat esnp_chr*.combined | awk '{OFS="\t"; $2=int($2/0.05); $3=($3>10)?10:$3; $4=($4>200000)?200000:$4; $4=int($4/10000);  print "MAF"$2"LD"$3"TSS"$4;}' | sort | uniq -c | awk '{OFS="\t"; print $2,$1}' > esnp.combined

# clean up
rm esnp_chr*

bsub -J "random[1-1000]" -q vshort -n 1 bash $pipeline_path/src/_deSNP.random.sh

## ---------------------------
## step 6: emperical p-value of N_deSNP = P(N_deSNP < n_deSNP)
## ---------------------------
Rscript $pipeline_path/modules/dSNP.eSNP.enrichment.R N_deSNP.txt `ls n_deSNP*.txt`  # output: dSNP.eSNP.enrichment.xls
rm n_deSNP*.txt