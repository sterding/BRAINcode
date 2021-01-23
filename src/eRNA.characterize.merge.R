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

# add class
df=subset(features, select=c(f06.TFBS, f07.P300, f08.CAGEenhancer, f09.chromHMM_brain, f12.DNaseROADMAP, f15.HCNE))
df$f06.TFBS=ifelse(df$f06.TFBS>=5,1,0)
df[df>0]=1;
features$class=ifelse(apply(df,1,sum)==0, 3, ifelse(df$f12.DNaseROADMAP==1, 1, 2))

write.table(features, "eRNA.characterize.xls", sep="\t", quote = F, col.names = NA, row.names = T)
