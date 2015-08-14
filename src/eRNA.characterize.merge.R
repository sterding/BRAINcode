# ================================
### merge features into one big file
## Usage: Rscript $HOME/neurogen/pipeline/RNAseq/src/eRNA.characterize.merge.R `ls eRNA.f*.txt`
# ================================

args<-commandArgs(TRUE)

## merge into one big file
n=length(args)

features=read.table(args[1], header=F);
rownames(features)=features[,1];

for(i in 2:n){
    message(paste0("[Merging file ", args[i], " ...] %", round(100*i/n, 1), " Done"));
    df=read.table(args[i], header=F); rownames(df) = df[,1];
    features=cbind(features, df[match(rownames(features), rownames(df)), 2])
}

features=features[,-1];
colnames(features)=gsub("eRNA.([^/]*).txt", "\\1", args[1:n])

write.table(features, "eRNA.characterize.xls", sep="\t", quote = F, col.names = NA, row.names = T)
