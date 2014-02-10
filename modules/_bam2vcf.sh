# Script to discover variations from input BAM file
# author: shuilin jin
# Date: 2014-1-24
# Version: 4.0
# Usage: bsub _bam2vcf.sh <input.bam>
a=/PHShome/sj750/neurogen/rnaseq_PD/run_output/HC_BN08-44_3
rna_seq_snp_calling=/PHShome/sj750/neurogen/sjn/rna_snp_calling/HC_BN08-44_3


# change MQ255_100

samtools view -h $a/accepted_hits.bam > $rna_seq_snp_calling/accepted_hits.sam

awk '{if ($5==255)$5==4; FQS="\t"; print}' $rna_seq_snp_calling/accepted_hits.sam > $rna_seq_snp_calling/accepted_hits_100.sam

samtools view -Sb $rna_seq_snp_calling/accepted_hits_100.sam -o $rna_seq_snp_calling/accepted_hits_100.bam


# add bam index:  (added after mark duplicate)

samtools index $rna_seq_snp_calling/accepted_hits_100.bam


# mark duplicate of bam file

java -Xmx50g -jar /PHShome/sj750/neurogen/local/picard/1.538/bin/MarkDuplicates.jar INPUT=$rna_seq_snp_calling/accepted_hits_100.bam OUTPUT=$rna_seq_snp_calling/accepted_hits_100_marked.bam METRICS_FILE=/PHShome/sj750/neurogen/sjn/rna_snp_calling/rna_seq_snp_calling/metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true


#realignerTargetCreator

java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T RealignerTargetCreator -U ALLOW_N_CIGAR_READS -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $rna_seq_snp_calling/accepted_hits_100_marked.bam -o $rna_seq_snp_calling/indel.intervals -known /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf 


# IndelRealigner

java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T IndelRealigner -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $rna_seq_snp_calling/accepted_hits_100_marked.bam -targetIntervals $rna_seq_snp_calling/indel.intervals -U ALLOW_N_CIGAR_READS -known /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o $rna_seq_snp_calling/indelrealigned_accpeted_hits_100_marked.bam



#add read group and sort

java -Xmx50g -jar /PHShome/sj750/projects/convert_bam/AddOrReplaceReadGroups.jar INPUT=$rna_seq_snp_calling/indelrealigned_accpeted_hits_100_marked.bam OUTPUT=$rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam RGPL=illumina RGPU=1 RGLB=bar RGSM=march2013111 RGID=foo


#index bam file
samtools index $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam


# HaplotypeCaller(after indelrealigner filter N Cigar)
java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T HaplotypeCaller -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --filter_reads_with_N_cigar -I $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam -o $rna_seq_snp_calling/SNPs_afterindelrealigner_filter_Ncigar.vcf -D /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/dbsnp_138.hg19.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -L /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v13.annotation.gtf.exons.bed


# HaplotypeCaller(after indelrealigner with N cigar
java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T HaplotypeCaller -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -U ALLOW_N_CIGAR_READS -I $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam -o $rna_seq_snp_calling/SNPs_afterindelrealigner_ALLOW_Ncigar.vcf -D /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/dbsnp_138.hg19.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -L /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v13.annotation.gtf.exons.bed


#BaseRecalibrator

java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T BaseRecalibrator -U ALLOW_N_CIGAR_READS -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam -knownSites /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/dbsnp_138.hg19.vcf -o $rna_seq_snp_calling/recal_data.table




#Print reads
java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T PrintReads -U ALLOW_N_CIGAR_READS -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked.bam -BQSR $rna_seq_snp_calling/recal_data.table -o $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked_recal.bam



# HaplotypeCaller(after BQSR)
java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T HaplotypeCaller -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --filter_reads_with_N_cigar -I $rna_seq_snp_calling/indelrealigned_accpeted_rg_hits_100_marked_recal.bam -o $rna_seq_snp_calling/SNPs_afterindelrealigner_afterBQSR.vcf -D /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/dbsnp_138.hg19.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -L /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v13.annotation.gtf.exons.bed




# Recalibration

java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T VariantRecalibrator -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --input $rna_seq_snp_calling/SNPs_afterindelrealigner_afterBQSR.vcf --maxGaussians 4 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /PHShome/sj750/neurogen/sjn/snp_calling/gatk_package/snp_calling/hg19/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $rna_seq_snp_calling/final/recalibrate.SNP.recal -tranchesFile $rna_seq_snp_calling/final/recalibrate.SNP.tranches -rscriptFile $rna_seq_snp_calling/final/recalibrate.SNP.plots.R



#Apply the desired level of recalibration to the SNPs in the call set

java -Xmx50g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T ApplyRecalibration -R /PHShome/sj750/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --input $rna_seq_snp_calling/SNPs_afterindelrealigner_afterBQSR.vcf -mode SNP --ts_filter_level 99.0 -recalFile $rna_seq_snp_calling/final/recalibrate.SNP.recal -tranchesFile $rna_seq_snp_calling/final/recalibrate.SNP.tranches -o $rna_seq_snp_calling/recalibrated_snps_afterBQSR_afterIndelrealigned_afterduplicate.vcf


