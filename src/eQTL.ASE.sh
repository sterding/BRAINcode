## ASE
cd ~/eRNAseq/HCILB_SNDA/
##################################################################
# step1: to generate VCF file using samtools mpileup for eQTL loci
##################################################################
#bam list
more +2 eQTL.covs.tab | awk '{print "/data/neurogen/rnaseq_PD/run_output/"$2"/accepted_hits.bam"}' > ASE.bamlist
#eQTL position
#cut -f1,3 /PHShome/xd010/eRNAseq/HCILB_SNDA/final.cis.eQTL.GWAS.FDR5pt.TFBS.xls  | sort -u > ASE.eQTLpos
cut -f1,3 final.cis.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.xls | sort -u > ASE.eQTLpos
bsub -q big-multi -n 4 -M 10000 -R 'rusage[mem=10000]' samtools mpileup -ABIs -q 0 -Q 0 -b ASE.bamlist -f $GENOME/Sequence/WholeGenomeFasta/genome.fa -l ASE.eQTLpos -o ASE.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.mpileup
# only include bi-allelic heterozygous SNPs (actually it's not necessary, as our SNPs already pass the MAF>5% filter, which means whatever left are heterozygous)
# vcftools --vcf SNP.vcf --min-alleles 2 --max-alleles 2 --out SNP.biallelic

##################################################################
# step2: run allelecounter to count # of reads for ALT/REF allele
# Note: https://github.com/secastel/allelecounter
##################################################################
SNP="/PHShome/gl871/neorogen/genotyping_PDBrainMap/eQTLMatrixBatch123/PDMAPHC.impute.vcf.gz"
module load pysam/0.8.4
module load htslib/1.2.1
#python ~/pipeline/modules/allelecounter.py --vcf $SNP --mpileup ASE.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.mpileup --sample "HC_UWA479" --samplelist ASE.bamlist
for i in `awk '$2==1' SNP.rs1664261 | cut -f1`; do python ~/pipeline/modules/allelecounter.py --vcf $SNP --mpileup ASE.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.mpileup --sample "$i" --samplelist ASE.bamlist; done | grep 61633127 > ASE.SNP.rs1664261

R
df=read.table("ASE.SNP.rs1664261", header=F, sep=" ")
p=wilcox.test(df$V7, df$V8, alternative = 'less', paired = T)$p.value
pdf("ASE.SNP.rs1664261.pdf", width=6, height=8)
barplot(height = -df$V7,add = F,axes = T, ylim=range(-df$V7, df$V8), names.arg=df$V1, las=2, cex.names=.7,col='blue', main=paste0("ASE on SNP rs1664261\n(paired wilcox p-value = ",round(p,2),")"), ylab='RNA-seq reads for heterozygous alleles')
barplot(height =  df$V8,add = T,axes = F, col='red')
legend('topleft',c('G allele','A allele'), col=c('red','blue'),pch=15)
dev.off()