#!/bin/bash

#Change headers to get FPKM
#Load modules
module load samtools-0.1.18

#Author: Robert Wang
#Rehead new file
#From @SQ	SN:1 to @SQ	SN:chr1
#Maps process

samtools view -H accepted_hits.bam | sed -E 's/SN:([0-9]+)\t/SN:chr\1\t/g;s/SN:([XY])\t/SN:chr\1\t/g;s/SN:MT\t/SN:chrM\t/g' > header.new

samtools reheader header.new accepted_hits.bam > new_bam
rm header.new
rm accepted_hits.bam
mv new_bam accepted_hits.bam

samtools view -h -o accepted_hits.sam accepted_hits.bam 
