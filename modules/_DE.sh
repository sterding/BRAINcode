###########################################
## bash script for running DE analysis using count data and covariance table
###########################################
#!/bin/sh

## HC vs. AD in TCPY
cd ~/neurogen/rnaseq_PD/results/merged

# begin R script
library("tidyverse")
library(RCurl)
setwd("~/neurogen/rnaseq_PD/results/merged")
sample_table_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vTNGgyiJHx9YItTT3bkI9F4VqYbpbrahX0RPQTiUBhtirwMhUc5bGoTYwbnkLlEkwkTkp1o9TEFQi3o/pub?gid=0&single=true&output=tsv"
subject_table_url="https://docs.google.com/spreadsheets/d/e/2PACX-1vTNGgyiJHx9YItTT3bkI9F4VqYbpbrahX0RPQTiUBhtirwMhUc5bGoTYwbnkLlEkwkTkp1o9TEFQi3o/pub?gid=16&single=true&output=tsv"
sampleTable=read.delim(textConnection(getURL(sample_table_url)), stringsAsFactors = F)
subjectTable=read.delim(textConnection(getURL(subject_table_url)), stringsAsFactors = F)
subjectTable = select(subjectTable, SOURCE_SUBJECT_ID, Source, Age,	Sex,	PMI, RIN, B.B,	CERAD,	NIA.R,	ABC.score)
filter(sampleTable, Cell_Type == 'TCPY') %>% 
  select(SAMPLE_ID=Sample_Name, SUBJECT_ID=Subject_ID, CONDITION=Diagnosis,	Batch) %>%
  mutate(Source=subjectTable$Source[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         Age = subjectTable$Age[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         Sex = subjectTable$Sex[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         PMI = subjectTable$PMI[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         RIN = subjectTable$RIN[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         #BB = subjectTable$B.B[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         #CERAD = subjectTable$CERAD[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         #NIA = subjectTable$NIA.R[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)],
         #ABC = subjectTable$ABC.score[match(SUBJECT_ID, subjectTable$SOURCE_SUBJECT_ID)]
         ) %>% 
  write.table(file="TCPY.covariates.txt", sep="\t", col.names = TRUE, row.names = F, quote = F)
# R end

Rscript ~/pipeline/modules/_DE.R -i genes.htseqcount.cufflinks.TCPY.uniq.xls -c TCPY.covariates.txt -o TCPY
