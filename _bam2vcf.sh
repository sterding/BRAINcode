# Script to discover variations from input BAM file
# author: shuilin jin
# Date: 2013-11-5
# Version: 1.0
# Usage: bsub _bam2vcf.sh <input.bam>

rna_seq_bamfile=${/PHShome/sj750/neurogen/sjn/hc_hg15/march2013192/accepted_hits.bam}
vcffile=${/PHShome/sj750/projects/snp_calling/march2013192/resultSNPs.vcf}

# Add RD group to RAN-seq bam file:

java -Xmx2g -jar /PHShome/sj750/projects/convert_bam/AddOrReplaceReadGroups.jar SORT_ORDER=coordinate INPUT=rna_seq_bamfile OUTPUT=/PHShome/sj750/projects/convert_bam/march2013192/accepted_hits.bam RGPL=illumina RGPU=1 RGLB=bar RGSM=march2013192 RGID=foo CREATE_INDEX=True



# sort the bam file by karyotypic 

java -Xmx5g -jar /PHShome/sj750/neurogen/local/picard/1.538/bin/ReorderSam.jar I=/PHShome/sj750/projects/convert_bam/march2013192/accepted_hits.bam O=/PHShome/sj750/projects/convert_bam/march2013192/kayrotypic.bam REFERENCE=/PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa



# add bam index:

samtools index /PHShome/sj750/projects/convert_bam/march2013192/kayrotypic.bam



# call snp

java -Xmx5g -jar /PHShome/sj750/projects/gatk_package/GenomeAnalysisTK.jar -T HaplotypeCaller -R /PHShome/sj750/neurogen/local/referenceGenome/hg19bt2/hg19.fa  --filter_reads_with_N_cigar -I /PHShome/sj750/projects/convert_bam/march2013192/kayrotypic.bam -o vcffile  -stand_call_conf 30.0 -stand_emit_conf 10.0









