## ASE
cd ~/eRNAseq/HCILB_SNDA/
querySNP=$1 # querySNP='rs1684902'

##################################################################
# step1: to generate VCF file using samtools mpileup for eQTL loci
##################################################################
#bam list
[ -e  ASE.bamlist ] || more +2 eQTL.covs.tab | awk '{print "/data/neurogen/rnaseq_PD/run_output/"$2"/accepted_hits.bam"}' > ASE.bamlist
## SNP position file
# #cut -f1,3 /PHShome/xd010/eRNAseq/HCILB_SNDA/final.cis.eQTL.GWAS.FDR5pt.TFBS.xls  | sort -u > ASE.eQTLpos
# cut -f1,3 final.cis.eQTL.d0.p1e-2.classI_or_nTFgt5.TFBS.xls | sort -u > ASE.eQTLpos

fgrep -w $querySNP final.cis.eQTL.d1e6.p1e-2.xls -m1 | awk '{OFS="\t"; print $8, $7;}' > ASE.$querySNP.pos

bsub -q big-multi -n 4 -M 10000 -R 'rusage[mem=10000]' samtools mpileup -ABIs -q 0 -Q 0 -b ASE.bamlist -f $GENOME/Sequence/WholeGenomeFasta/genome.fa -l ASE.$querySNP.pos -o ASE.$querySNP.mpileup

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
head -n1 ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt > temp
fgrep -m1 $querySNP ~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt >> temp  ## only the first hit
rowsToCols temp ASE.$querySNP.subject
for i in `awk '$2==1' ASE.$querySNP.subject | cut -f1`; do python ~/pipeline/modules/allelecounter.py --vcf $SNP --mpileup ASE.$querySNP.mpileup --sample "$i" --samplelist ASE.bamlist; done > ASE.SNP.$querySNP

R
querySNP='rs1684902'
df=read.table(paste("ASE.SNP",querySNP,sep="."), header=F, sep=" ")
p=wilcox.test(df$V7, df$V8, alternative = 'less', paired = T)$p.value
pdf(paste("ASE.SNP",querySNP,"pdf",sep="."), width=6, height=8)
barplot(height = -df$V7,add = F,axes = T, ylim=range(-df$V7, df$V8), names.arg=df$V1, las=2, cex.names=.7,col='blue', main=paste0("ASE on SNP ", querySNP,"\n(paired wilcox p-value = ",round(p,2),")"), ylab='RNA-seq reads for heterozygous alleles')
barplot(height =  df$V8,add = T,axes = F, col='red')
legend('topleft',c('G allele','A allele'), col=c('red','blue'),pch=15)

# percentage plot
barplot(colMeans(df[,7:8]>0), ylim=c(0,.5), col=c('red','blue'), names.arg = c('A','G'))

dev.off()