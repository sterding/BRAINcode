# Script to discover variations from input BAM file
# author: shuilin jin
# Date: 2013-11-5
# Version: 1.0
# Usage: bsub _bam2vcf.sh <input.bam>

rna_seq_bamfile=$1
vcffile=${/name/bam$/vcf}

# Add RD group to RNA-seq bam file:

java -Xmx5g -jar /PHShome/sj750/projects/convert_bam/AddOrReplaceReadGroups.jar SORT_ORDER=coordinate INPUT=$rna_seq_bamfile OUTPUT=rna_seq_bamfile_ad.bam RGPL=illumina RGPU=1 RGLB=bar RGSM=march2013192 RGID=foo CREATE_INDEX=True



# sort the bam file by karyotypic 

java -Xmx5g -jar /PHShome/sj750/neurogen/local/picard/1.538/bin/ReorderSam.jar I=rna_seq_bamfile_ad.bam O=rna_seq_bamfile_ad.kayrotypic.bam REFERENCE=/PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa



# add bam index:

samtools index rna_seq_bamfile_ad.kayrotypic.bam




## mark duplicate of bam file
java -Xmx5g -jar /PHShome/sj750/neurogen/local/picard/1.538/bin/MarkDuplicates.jar I=/PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march2013101_first.bam O=/PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march2013101_marked.bam METRICS_FILE=/PHShome/sj750/neurogen/ranSeq_snpcalling/march101/ metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true


##Analyze patterns of covariation in the sequence dataset

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T BaseRecalibrator -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa  -I /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_realigned_reads.bam -knownSites dbsnp.vcf   -knownSites gold_indels.vcf   -o /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_data.table 

##Do a second pass to analyze covariation remaining after recalibration
#Action

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T BaseRecalibrator -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa -I /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_realigned_reads.bam -knownSites dbsnp.vcf   -knownSites gold_indels.vcf -BQSR /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_data.table   
 -o /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_post_recal_data.table

##Generate before/after plots

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T AnalyzeCovariates  -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa
-before /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_data.table    -after /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_post_recal_data.table -plots /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recalibration_plots.pdf


##Apply the recalibration to your sequence data

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar 
    -T PrintReads 
    -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa
    -I /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_realigned_reads.bam
    -BQSR /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_data.table    -o /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_reads.bam

## Compress read data with ReduceReads

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T ReduceReads   -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa  -I /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_recal_reads.bam      -o /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_reduced_reads.bam 


## merge all the bam file together….

samtools merge out.bam in1.bam in2.bam in3.bam



## Call variants on a diploid genome with the HaplotypeCaller 

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar
    -T HaplotypeCaller 
    -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa
    -I /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_reduced_reads.bam 
    --genotyping_mode DISCOVERY 
    -stand_emit_conf 10 
    -stand_call_conf 30  
    -o /PHShome/sj750/neurogen/ranSeq_snpcalling/march101/march101_raw_variants.vcf 



# call snp

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T HaplotypeCaller -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa  --filter_reads_with_N_cigar -I rna_seq_bamfile_ad.kayrotypic.bam -o vcffile  -stand_call_conf 30.0 -stand_emit_conf 10.0









