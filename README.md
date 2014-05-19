RNA-seq Analysis Pipeline
=========================

File convention
---------------
To standardlize the process, we suggest to start with conventional file name format. In this pipeline, we will use filename convention as:
`[PD|AD|HD|HC]_[brainBankID]_<batchID>_<rep#>.[R1|R2].<fq|fastq><.gz>`
where,
- Beginning with patient type, PD="Parkinson's disease", AD="Alzhimer's disease", HD="Huntington's disease" and HC="Health controls", in 2-letter abbreviation;
- Following with unique BrainBankID rather than the lane# or date
- batchID and rep# are optional, in format of (if any) 
 - batchID: batch1, batch2 etc.
 - rep#: rep1 is for technique replicate, Rep1 for biological replicate
- Use R1 or R2 to tell the two mates of pair-end sequencing data
- Use unzipped or zipped fastq, but NOT both.

Folder structure
----------------
For each project (e.g. PDMap), it should have a folder named with project short name and followed by the type of data, such as PDMap_RNAseq, HDPredict_smallRNA etc. 
- rawfiles
 - Raw files from sequencing;
 - Making soft links for conventional file name
 - Log file (readme.txt or Excel)
- filtered
 - Filtered files (e.g. adaptor removal/clip)
- run_output 
 - sample1
  - Output of Tophat/Cufflinks/htseq-count runs 
 - sample2 etc. 
- for_display
 - Files used for display on UCSC / IGV, such as *.bam, *.bam.bai, *.bw, *.gtf etc.
 - Can be soft links to the output files

Pipeline Requirement
--------------------
1. Install programs: tophat, cufflinks, bowtie, bedtools, samtools, fastq-mcf, fastqc, and htseq-count; and add their executable programs in the $PATH

Pipeline structure
--------------
### Main script:
RNAseq.pipeline.sh
### Modules:
#### _RNAseq.sh

#### _bam2vcf.sh
#### _callSNP.sh
#### _sam2variation.awk

#### _clustComRNASeq.R
performs sample clustering analysis based on genes.fpkm_tracking and isoforms.fpkm_tracking from Cufflinks' output  
#### _cluster.sh

#### _DE_cuffdiff.sh
#### _DE_DEseq.R

#### _factor_analysis.R

