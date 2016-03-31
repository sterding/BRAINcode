# Script to conduct enrichment analysis of disease SNPs and eQTL SNPs

#Usage: 

eQTL_final=$1
eQTL_final=~/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/final.cis.eQTL.xls

dir="${eQTL_final%/*}"

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

## ---------------------------
## step 3: calculate number of eSNPs that are also GWAS dSNPs (N_deSNP)
## ---------------------------
cd $dir
SNPpos=/data/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID
# p-value < 0.001
awk '$5<=0.001' $eQTL_final | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' > ${eQTL_final/xls/bed}
## FDR < -.05
#awk '$6<=0.05' $eQTL_final | cut -f1 | sort -u | fgrep -f - $SNPpos | awk '{OFS="\t"; print $2,$3-1,$3,$1,$4}' > ${eQTL_final/xls/bed}
intersectBed -a ${eQTL_final/xls/bed} -b $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pruned.bed -wo | cut -f12 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > N_deSNP.txt

## ---------------------------
## step 4: generate 1000 random sets of SNPs with matched MAF, number of LD partners (r2>0.8), and distance from transcription start site (TSS) with eSNP. For each set, calculate number of each eSNP set that are also dSNPs (n_deSNP)
## ---------------------------
# disect eSNP into grids by MAF, number of LD partners (r2>0.8), and distance from transcription start site (TSS)
GENEpos=$GENOME/Annotation/Genes/gencode.v19.annotation.gtf.genes.bed
for i in `seq 1 22`; do
  echo $i;
  # only SNPs on the array
  awk -vchr=$i '{OFS="\t"; if($2=="chr"chr) print $2,$3-1,$3,$1}' $SNPpos | intersectBed -a <(awk '{OFS="\t"; print "chr"$1,$4-1,$4,$2}' $GENOME/Annotation/Variation/1000G/Chr$i.bim) -b - -u > esnp_chr$i.bed
  cut -f4 esnp_chr$i.bed > esnp_chr$i
  bsub -q short -n 1 plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --freq --extract esnp_chr$i --silent --out esnp_chr$i
  bsub -q normal -n 1 plink --bfile $GENOME/Annotation/Variation/1000G/Chr$i --show-tags esnp_chr$i --list-all --tag-kb 250 --tag-r2 0.8 --silent --out esnp_chr$i
  bedtools closest -a esnp_chr$i.bed -b <(sortBed -i $GENEpos) -d -t first | awk '$11==0' | awk '{OFS="\t"; tss=($10=="+")?$6:$7; d=tss-$2; if(d>0) d=-d; print $4, d;}' > esnp_chr$i.d2TSS
  bedtools closest -a esnp_chr$i.bed -b <(sortBed -i $GENEpos) -d -t first | awk '$11!=0' | cut -f1-4 | closestBed -a - -b <(awk '{OFS="\t"; tss=($6=="+")?$2:($3-1);  print $1, tss, tss+1, $4, $3-$2, $6}' $GENEpos | sortBed) -D b -t first | awk '{OFS="\t"; print $4,($11<0)?-$11:$11;}' >> esnp_chr$i.d2TSS
done
# combine
for i in `seq 1 22`; do paste <(awk '{OFS="\t"; if(NR>1) print $2,$5}' esnp_chr$i.frq | sort -k1,1) <(awk '{OFS="\t"; if(NR>1) print $1,$4}' esnp_chr$i.tags.list | sort -k1,1 | cut -f2) <(sort -k1,1 esnp_chr$i.d2TSS | cut -f2) > esnp_chr$i.combined; done
# direct
cat esnp_chr*.combined | awk '{OFS="\t"; $2=int($2*10); $3=($3<=3)?1:{($3<=10)?2:3}; $4=int($4/10000); print;}' > esnp.combined

## step 5: emperical p-value of N_deSNP = P(N_deSNP < n_deSNP)


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
echo "# eRNA-private"
[ -e eRNA.private.major.bed ] && intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.private.major.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.private.major.HTNE
[ -e eRNA.private.minor.bed ] && intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.private.minor.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.private.minor.HTNE
echo "# mRNA inner exons"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $GENOME/Annotation/Genes/mRNA.innner.exon.bed  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.exon
echo "# promoter -- [-200,+200] of protein-coding GENCODE v19 TSS"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.promoter
echo "# randomly sampling"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.random.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.random

#rm $snps_in_LD.autosomal.associations.bed

echo "## total overlapped SNPs count with each dataset"
### ##################
echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.$type.counts.summary
echo "HTNE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "exon" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $GENOME/Annotation/Genes/mRNA.innner.exon.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "promoter" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b $GENOME/Annotation/Genes/gencode.v19.annotation.pc.promoter.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "random" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.random.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary

echo "## Fisher test and make plot"
### ##################
Rscript $pipeline_path/src/eRNA.SNP.enrichment.R $type