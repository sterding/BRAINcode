# Usage: 
# for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh PLINK $i; done 
# for i in HCILB_SNDA HC_PY HC_nonNeuron; do echo $i; bsub -n 1 -q normal -J $i bash $pipeline_path/src/eRNA.SNP.enrichment.sh SNAP $i; done 

type=$1
samplegroup=$2

cd ~/eRNAseq/$samplegroup

## pre-steps to get LD data for GWAS SNPs: 
## Ref ~/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/GWASCatalog/README.txt

[ "$type" == "SNAP" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.SNAP.LD_w250.r2_0.8.bed
[ "$type" == "PLINK" ] && snps_in_LD=$GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.snps_in_LD.PLINK.LD_w250.r2_0.8.bed

## extract all autosomal.associations
awk 'BEGIN{FS="\t"; OFS="\t";}{split($8,a,"|");  n=split(a[2],b,";"); print $1,$2,$3,$7; for(i=1;i<n;i++) print $1,b[i]-1,b[i],$7;}' $snps_in_LD | grep -v chrX | grep -v chrY | sortBed | uniq > $snps_in_LD.autosomal.associations.bed

# number of gwas SNPs
wc -l $snps_in_LD
# number of associations
wc -l $snps_in_LD.autosomal.associations.bed

# random regions
#bedtools shuffle -excl <(cat blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -noOverlapping -i eRNA.bed -g $GENOME/Annotation/Genes/ChromInfo.txt > eRNA.random.bed
bedtools random -n 100000 -l 400 -seed 1234 -g $GENOME/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(cat ../blacklist.bed eRNA.bed | cut -f1-3 | sortBed | mergeBed -i -) -v > eRNA.random.bed
# wc -l eRNA.random.bed --> N=72202 (close to 71469 lines of eRNA.bed)

echo "## overlapped SNPs with each dataset"
### ##################
echo "# all"
cat $snps_in_LD.autosomal.associations.bed | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.all
echo "# eRNA"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.HTNE
echo "# eRNA-private"
[ -e eRNA.private.bed ] && intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.private.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.pHTNE
echo "# mRNA inner exons"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b <(grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed | awk '{if(id!=$4) id=$4; else print}' | sort -k4,4 -k2,2nr | awk '{if(id!=$4) id=$4; else print}')  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.exon
echo "# promoter -- [-200,+200] of protein-coding GENCODE v19 TSS"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b <(grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400}')  -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.promoter
echo "# randomly sampling"
intersectBed -a $snps_in_LD.autosomal.associations.bed -b eRNA.random.bed -u | cut -f4 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > SNP.$type.count.random

#rm $snps_in_LD.autosomal.associations.bed

echo "## total overlapped SNPs count with each dataset"
### ##################
echo "all" `wc -l $GENOME/Annotation/Variation/snp137.bed.groupped.SNP | cut -f1 -d' '` > SNP.$type.counts.summary
echo "HTNE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
[ -e eRNA.private.bed ] && echo "pHTNE" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.private.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "exon" `grep protein_coding.protein_coding $GENOME/Annotation/Genes/exons.bed | awk '{if(id!=$4) id=$4; else print}' | sort -k4,4 -k2,2nr | awk '{if(id!=$4) id=$4; else print}' | intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b stdin -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "promoter" `grep protein_coding.protein_coding $GENOME/Annotation/Genes/gencode.v19.annotation.bed12 | awk '{OFS="\t"; s=($6=="+")?($2-200):($3-200); if(s<0) s=0; print $1,s,s+400}'| intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b stdin -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary
echo "random" `intersectBed -a $GENOME/Annotation/Variation/snp137.bed.groupped.SNP -b eRNA.random.bed -u | wc -l | cut -f1 -d' '` >> SNP.$type.counts.summary

echo "## Fisher test and make plot"
### ##################
Rscript $pipeline_path/src/eRNA.SNP.enrichment.R $type