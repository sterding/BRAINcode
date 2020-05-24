##############################################
## R script to visualize interaction data
## Author: Xianjun Dong
## Usage: Rscript $0 assocation.txt annotation.txt
## Date:2017-05-12
## Version: 0.0 (Unfinished)
############################################

args<-commandArgs(TRUE)
fileAssociation = args[1]  
fileAnnotation  = args[2]
type = args[3]

# Format of assocaition file:
# x1 x2 0.9
# x1 x3 0.8
# x2 x4 0.3
# ...

# Format of annotation file:
# chr1  123 456 x1
# chr1  789 496 x2
# chr1  013 856 x3
# ...

# other optional parameters
type='eQTL' # or 'HiC', 'LD', 'coexpression'

#################
## 1. read data
#################
df  <- read.table(fileAssociation, header=F, col.names = c("name1","name2","value"))
pos <- read.table(fileAnnotation, header=F, col.names = c("chr","start","end","name"))

#################
## 2. proceed the data
#################

