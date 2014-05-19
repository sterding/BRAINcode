#################################
# coverage of genome by RNAseq 
#################################
# % of genome covered by RNAseq (>5rpm)
DONE
# Rscript to make plot of coverage vs. rpm
RPM=data.frame();
pdf("rpm_vs_coverage.merged.pdf", width=6, height=6)
for(i in list.files(pattern="*rpm_vs_coverage.txt")){
    samplename=gsub(".accepted_hits", "", gsub(".normalized.rpm_vs_coverage.txt","",i))
    df=read.table(i, header=F); rownames(df)=df[,1]; df=df[,-1,drop=F]; colnames(df)=samplename
    df[,1]=cumsum(df[,1])*100/3e9
    if(ncol(RPM)==0) {RPM=df;} else {RPM=cbind(RPM, df);}
    plot(df[,1], xlab="RPM", ylab="% of cummulative base pairs covered", xaxt='n',yaxt='n', type="l", main=samplename, log="y", ylim=c(.005,100))
    axis(side=1, at=c(0:10)*10+1, labels=rownames(RPM)[c(0:10)*10+1])
    axis(side = 2, at = (locs <- 100/c(1,10,100,1000)), labels = locs)
    axis(side = 2, at = (locs2 <- c(outer(1:10, c(0.1, 1, 10, 100), "/"))), labels = NA, tcl = -0.3)
}
# merge of all
plot(RPM[,1], xlab="RPM", ylab="% of cummulative base pairs covered", xaxt='n',yaxt='n', main="all together", type="n", log="y", ylim=c(.005,100))
#sapply(1:ncol(RPM), function(x) lines(RPM[,x], col=x))
sapply(colnames(RPM), function(x) lines(RPM[,x], col=as.integer(sub(".*_([1-4]).*","\\1", x))))
points(rep(nrow(RPM),ncol(RPM)), RPM[nrow(RPM),], col=as.integer(sub(".*_([1-4]).*","\\1", colnames(RPM))))
axis(side=1, at=c(0:10)*10+1, labels=rownames(RPM)[c(0:10)*10+1])
axis(side = 2, at = (locs <- 100/c(1,10,100,1000)), labels = locs)
axis(side = 2, at = (locs2 <- c(outer(1:10, c(0.1, 1, 10, 100), "/"))), labels = NA, tcl = -0.3)
legend("topleft", paste("batch", 1:4), col=1:4, lty=1, cex=.8, bty='n')

dev.off()

echo -n "rpm" > ../results/coverage/rpm_vs_coverage.allsamples.txt;
for i in */accepted_hits.normalized.rpm_vs_coverage.txt; do echo -ne "\t"${i/\/*/} >> ../results/coverage/rpm_vs_coverage.allsamples.txt; done
echo '' >> ../results/coverage/rpm_vs_coverage.allsamples.txt
paste */accepted_hits.normalized.rpm_vs_coverage.txt | grep -v "#" | awk '{printf("%s", $1);i=2; while(i<=NF) {printf("\t%s",$i);i=i+2;} printf("\n");}' >> ../results/coverage/rpm_vs_coverage.allsamples.txt

RPM=read.table("/PHShome/xd010/neurogen/rnaseq_PD/results/coverage/rpm_vs_coverage.allsamples.txt", header=T);
rownames(RPM)=RPM[,1]; RPM=RPM[,-1, drop=F]

rpm=sapply(colnames(RPM), function(x) cumsum(RPM[,x])*100/3e9)
rownames(rpm)=rownames(RPM)
> summary(rpm["0.05",])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.9494  5.9520  6.9810  6.8110  8.2490 10.1500 
> summary(rpm["0",])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.706  12.840  17.330  16.610  20.510  30.210 

pdf("~/neurogen/rnaseq_PD/results/coverage/rpm_vs_coverage.pdf", width=6, height=6)
plot(RPM[,1], xlab="RPM", ylab="% of cummulative base pairs covered", xaxt='n',yaxt='n', main="all together", type="n", log="y", ylim=c(.005,100))
#sapply(1:ncol(RPM), function(x) lines(cumsum(RPM[,x])*100/3e9, col=x))
#legend("topleft", paste(colnames(RPM), "(",format(apply(RPM, 2, sum)*100/3e9,digits=3), ")"), col=1:ncol(RPM), lty=1, cex=.5, bty='n')
sapply(colnames(RPM), function(x) lines(cumsum(RPM[,x])*100/3e9, col=as.integer(sub(".*_([1-4])","\\1", x))))
legend("topleft", paste(colnames(RPM), "(",format(apply(RPM, 2, sum)*100/3e9,digits=3), ")"), col=1:ncol(RPM), lty=1, cex=.5, bty='n')
points(rep(nrow(RPM),ncol(RPM)), apply(RPM, 2, sum)*100/3e9, col=1:ncol(RPM))
axis(side=1, at=c(0:10)*10+1, labels=rownames(RPM)[c(0:10)*10+1])
axis(side = 2, at = (locs <- 100/c(1,10,100,1000)), labels = locs)
axis(side = 2, at = (locs2 <- c(outer(1:10, c(0.1, 1, 10, 100), "/"))), labels = NA, tcl = -0.3)
dev.off()

#################################
# combined bigwig, using unionBedGraphs in bedtools
#################################
mkdir ~/neurogen/rnaseq_PD/results/merged
cd ~/neurogen/rnaseq_PD/results/merged
> .paraFile;
# save the following code block to paraFile

## version 1: mean of all samples

# multi-mapper
#unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_MCPY.mean.multi.normalized.bedGraph && bedGraphToBigWig HC_MCPY.mean.multi.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_MCPY.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > PD_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig PD_SNDA.mean.multi.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt PD_SNDA.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig HC_SNDA.mean.multi.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_SNDA.mean.multi.normalized.bw
#unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > ILB_SNDA.mean.multi.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.mean.multi.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt ILB_SNDA.mean.multi.normalized.bw

# uniq-mapper
unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_MCPY.mean.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.mean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_MCPY.mean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > PD_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.mean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt PD_SNDA.mean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > HC_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.mean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_SNDA.mean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/*normalized.bedGraph` | awk '{OFS="\t"; s=0; for(i=4;i<=NF;i++) s=s+$i; s=s/(NF-3); print $1,$2,$3,s}' > ILB_SNDA.mean.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.mean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt ILB_SNDA.mean.uniq.normalized.bw

# version 2: median of all samples

unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, median(S)}' > HC_MCPY.median.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.median.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_MCPY.median.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, median(S)}' > PD_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.median.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt PD_SNDA.median.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, median(S)}' > HC_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.median.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_SNDA.median.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function median(v) {c=asort(v,j);  if (c % 2) return j[(c+1)/2]; else return (j[c/2+1]+j[c/2])/2.0; } {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, median(S)}' > ILB_SNDA.median.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.median.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt ILB_SNDA.median.uniq.normalized.bw

# version 3: trimmed mean (10%)

unionBedGraphs -i `ls ../../run_output/HC*_MCPY_[234]/uniq/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, trimmedMean(S, 0.1)}' > HC_MCPY.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig HC_MCPY.trimmedmean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_MCPY.trimmedmean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/PD*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, trimmedMean(S, 0.1)}' > PD_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig PD_SNDA.trimmedmean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt PD_SNDA.trimmedmean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/HC*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, trimmedMean(S, 0.1)}' > HC_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig HC_SNDA.trimmedmean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt HC_SNDA.trimmedmean.uniq.normalized.bw
unionBedGraphs -i `ls ../../run_output/ILB*_SNDA_[234]/uniq/*normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); s=0; for(i=a+1;i<=(c-a);i++) s+=j[i];return s/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; print $1,$2,$3, trimmedMean(S, 0.1)}' > ILB_SNDA.trimmedmean.uniq.normalized.bedGraph && bedGraphToBigWig ILB_SNDA.trimmedmean.uniq.normalized.bedGraph /data/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt ILB_SNDA.trimmedmean.uniq.normalized.bw

ParaFly -c .paraFile -CPU 8

#################################
# calculate rpm vs. coverage for individual sample
#################################
cd ~/neurogen/rnaseq_PD/results/coverage

ls ../../run_output/ | rowsToCols stdin all.uniq.accepted_hits.normalized.rpm_vs_coverage.txt
> tmp; for i in ../../run_output/*/uniq/accepted_hits.normalized.rpm_vs_coverage.txt; do echo $i; sed 's/#.*=/\t/' $i | cut -f2 | paste - tmp> tmp2; mv tmp2 tmp; done
cat tmp >> all.uniq.accepted_hits.normalized.rpm_vs_coverage.txt

> merge.pt05.bed
>merged.cov.txt
for i in ../../run_output/*/uniq/accepted_hits.normalized.*ph; do awk '$4>=0.05' $i | sed 's/ /\t/g' | cut -f1-3 | cat - merge.pt05.bed | sortBed | mergeBed > $i.merged.bed; cp $i.merged.bed merge.pt05.bed; echo $i; echo $i `awk 'BEGIN{s=0}{s=s+($3-$2);}END{ print s*100/3e9}' merge.pt05.bed` >> merged.cov.txt; done

> merge.5reads.bed
>merged.cov.5reads.txt
for i in ../../run_output/*/uniq/accepted_hits.bedGraph; do awk '$4>=5' $i | cut -f1-3 | cat - merge.5reads.bed | sortBed | mergeBed > $i.merged.bed; ln -f $i.merged.bed merge.5reads.bed; echo $i; echo $i `awk 'BEGIN{s=0}{s=s+($3-$2);}END{ print s}' merge.5reads.bed` >> merged.cov.5reads.txt; done

rm merge.pt05.bed merge.5reads.bed

## Rscript
pdf("merged.cov.rpm0.05.pdf", width=7,height=5)
par(mar=c(6,4,2,2))
df2=read.table("merged.cov.5reads.txt");
df=read.table("merged.cov.txt")
plot(df$V2,pch=16,cex=.6, xlab="", ylab="Percentage(%) of coverage by union of samples", xaxt = "n", ylim=range(df$V2, df2$V2*100/3e9));
points(df2$V2*100/3e9, pch=1, cex=.6)
axis(1, 1:nrow(df), gsub(".*run_output/(.*)/uniq.*", "\\1", df$V1), las=2, cex.axis=.5)
legend("topleft", c("RPM>=0.05", "at least 5 reads"),pch=c(16,1))
dev.off()