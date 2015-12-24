## ASE

##################################################################
# step1: to generate VCF file using samtools mpileup for eQTL loci
##################################################################
#bam list
more +2 ~/eRNAseq/HCILB_SNDA/eQTL.covs.tab | awk '{print "/data/neurogen/rnaseq_PD/run_output/"$2"/uniq/accepted_hits.bam"}' > ASE.bamlist
#eQTL position
cut -f1,3 /PHShome/xd010/eRNAseq/HCILB_SNDA/final.cis.eQTL.GWAS.FDR5pt.TFBS.xls > ASE.eQTLpos
bsub -q big-multi -n 4 -M 10000 -R 'rusage[mem=10000]' samtools mpileup -ABI -q 0 -Q 0 -b ASE.bamlist -f $GENOME/Sequence/WholeGenomeFasta/genome.fa -l ASE.eQTLpos -o ASE.eQTL.mpileup
# only include bi-allelic heterozygous SNPs
vcftools --vcf SNP.vcf --min-alleles 2 --max-alleles 2 --out SNP.biallelic

##################################################################
# step2: run allelecounter to count # of reads for ALT/REF allele
# Note: https://github.com/secastel/allelecounter
##################################################################

allelecounter --vcf SNP.biallelic.vcf --mpileup ASE.eQTL.mpileup --sample xxxx --o ASE.xxxx.out