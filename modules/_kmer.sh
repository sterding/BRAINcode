#!/bin/bash

sample=$1
K=$2

kpal count -p $sample -k 9 ${sample}.R1.fa ${sample}.R1.k9
kpal count -p $sample -k 9 ${sample}.R2.fa ${sample}.R2.k9
kpal merge ${sample}.R1.k9 ${sample}.R2.k9 ${sample}.k9