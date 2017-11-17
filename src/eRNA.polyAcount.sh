# script to check if the TNE is biasly selected by polyA enrichment (based on discussion with James Gusella)

inputBed=$1 # inputBed=eRNA.bed

bedtools slop -i $inputBed -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size -b 100 | bedtools getfasta -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -tab -bed - | awk '{OFS="\t"; seq=tolower($2); print $1,length(seq), gsub("a","a",seq), gsub("aa","aa",seq), gsub("aaa","aaa",seq),gsub("aaaa","aaaa",seq),gsub("aaaaa","aaaaa",seq),gsub("aaaaaa","aaaaaa",seq),gsub("aaaaaaa","aaaaaaa",seq),gsub("aaaaaaaa","aaaaaaaa",seq),gsub("aaaaaaaaa","aaaaaaaaa",seq),gsub("aaaaaaaaaa","aaaaaaaaaa",seq);}' > $inputBed.polyAcount

bedtools shuffle -seed 123 -noOverlapping -i $inputBed -g $GENOME/Annotation/Genes/ChromInfo.txt | bedtools slop -i - -g $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size -b 100 | bedtools getfasta -fi $GENOME/Sequence/WholeGenomeFasta/genome.fa -tab -bed - | awk '{OFS="\t"; seq=tolower($2); print $1,length(seq), gsub("a","a",seq), gsub("aa","aa",seq), gsub("aaa","aaa",seq),gsub("aaaa","aaaa",seq),gsub("aaaaa","aaaaa",seq),gsub("aaaaaa","aaaaaa",seq),gsub("aaaaaaa","aaaaaaa",seq),gsub("aaaaaaaa","aaaaaaaa",seq),gsub("aaaaaaaaa","aaaaaaaaa",seq),gsub("aaaaaaaaaa","aaaaaaaaaa",seq);}' > $inputBed.polyAcount_random

# R script
R
df1=read.table('eRNA.bed.polyAcount', header = F, sep = "\t", stringsAsFactors = F, row.names = 1, col.names = c("ID","length",strrep('A',1:10)))
df2=read.table('eRNA.bed.polyAcount_random', header = F, sep = "\t", stringsAsFactors = F, row.names = 1, col.names = c("ID","length",strrep('A',1:10)))

#sapply(1:11, function(i) t.test(df1[,i], df2[,i], paired = T)$p.value)

df=data.frame(TNE=colSums(df1), random=colSums(df2))
nuc_length = df[1,1]
df=df[-1,]
## convert to %
#df_P = t(as.matrix(df / (nuc_length- nrow(df1)*(1:10-1))))
#barplot(df_P, beside = TRUE,legend.text = rownames(df_P), args.legend = list(x = "topright", bty="n"),las=2, log='y')

pdf("eRNA.bed.polyAcount.pdf", width=5.5, height = 6)
par(mar=c(8,4,2,2));barplot(log10(t(as.matrix(df))),beside = TRUE,legend.text = rownames(df_P), args.legend = list(x = "topright", bty="n"), las=2,ylim = c(3.8,7.5),xpd = FALSE, yaxt='n', ylab="Number of polyA occurances", main=paste("P-value:",signif(t.test(df[,1],df[,2],paired = T)$p.value, 2)))
#minor.ticks.axis(2,9,mn=3.8,mx=7.5)
log10.axis <- function(side, at, ...) {
  at.minor <- log10(outer(1:9, 10^(min(at):(max(at)-1))))
  lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
  axis(side=side, at=at, labels=lab, ...)
}
log10.axis(2, at=seq(4, 7, 1))
dev.off()


