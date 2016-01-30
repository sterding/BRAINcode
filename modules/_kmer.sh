#!/bin/bash

sample=$1
K=9

kpal count -p $sample -k $K ${sample}.R1.fa ${sample}.R1.k$K
kpal count -p $sample -k $K ${sample}.R2.fa ${sample}.R2.k$K
kpal merge ${sample}.R1.k$K ${sample}.R2.k$K ${sample}.k$K
rm ${sample}.R1.k$K ${sample}.R2.k$K