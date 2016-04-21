#!/bin/bash
# usage: bsub -J "random[1-1000]" -q vshort -n 1 bash $pipeline_path/src/_deSNP.random.sh \$LSB_JOBINDEX
# either use $LSB_JOBINDEX in the script (ref: http://stackoverflow.com/questions/11212923/referencing-job-index-in-lsf-job-array)
# or, \$LSB_JOBINDEX outside of script as an argument (ref: https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_admin/job_array_cl_args.dita)

i=$LSB_JOBINDEX  
> random$i.list
while read pattern count; do shuf -n $count /data/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/$pattern >> random$i.list; done < esnp.combined
# get bed for random selection
cat /data/neurogen/rnaseq_PD/results/eQTL/HCILB_SNDA/ALL_chr*.bed | fgrep -w -f random$i.list | intersectBed -a - -b $GENOME/Annotation/GWASCatalog/gwas_catalog_v1.0-downloaded.hg19.pruned.bed -wo | cut -f11 | sed 's/ (.*//g;s/;.*//g;s/ /_/g' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2, $1}' > n_deSNP$i.txt
rm random$i.list
