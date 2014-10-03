# script to divide mapped reads into different annotation region
# usage: reads_distribution.sh input.bam
inputBAMs=$1  # wild card of bam files e.g. ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.bam
export ANNOTATION=$GENOME/Annotation/Genes

ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/exons.bed {}`' | sort | cut -f2 > /tmp/exons
ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/introns.bed {}`' | sort | cut -f2 > /tmp/introns
ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/5utr.bed {}`' | sort | cut -f2 > /tmp/5utr
ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/3utr.bed {}`' | sort | cut -f2 > /tmp/3utr
ls $inputBAMs | parallel 'echo {} `samtools view -c {} chrM`' | sort | cut -f2 > /tmp/chrM
ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/LINE.bed {}`' | sort | cut -f2 > /tmp/LINE
ls $inputBAMs | parallel 'echo {} `samtools view -c -L $ANNOTATION/SINE.bed {}`' | sort | cut -f2 > /tmp/SINE

ls $inputBAMs | paste -d' ' - /tmp/5utr /tmp/exons /tmp/introns /tmp/3utr /tmp/chrM /tmp/LINE /tmp/SINE 

