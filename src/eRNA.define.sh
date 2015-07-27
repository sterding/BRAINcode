## script to define HiTNE with the following features
# Usage: $pipeline_path/eRNA.define.sh HCILB_SNDA
# Author: Xianjun Dong
# Date: July 17, 2015
# preqisition: 
# 1. run bash /data/neurogen/pipeline/RNAseq/src/get.regulatoryMarks.data.sh to Download bigwig files for various marks
# 2. run bash /data/neurogen/pipeline/RNAseq/src/get.background.sh to generate the background region

## TODO: modify code to work for any sample group

################################################
# HiTNE definition:
# 1) density higher than the basal level,  
# 2) summit >0.05 RPM, --> p<0.05 comparing to the transcriptional noise
# 3) located in non-generic regions (e.g. 500bp away from any annotated exons),
# 4) at least 100bp in length,
# 5) don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly) 
# 6) q-value<0.05 in at least 25% samples when comparing with random non-functional background
################################################

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

cd ~/projects/PD/results/eRNA

SAMPLE_GROUP=$1
inputBG=/data/neurogen/rnaseq_PD/results2/merged/trimmedmean.uniq.normalized.$SAMPLE_GROUP.bedGraph

mkdir $SAMPLE_GROUP
cd $SAMPLE_GROUP

# ===========================================================================
#: background region to measure transcriptional noise: genomic regions excluding the known regions with RNA activities (known exons+/-500bp, rRNA, CAGE-defined enhancers, promoters)
# ===========================================================================

# RNAseq signal distribution in the background region
intersectBed -a $inputBG -b ../blacklist.bed -sorted -v | awk '{OFS="\t"; print $3-$2, $4}' | shuf -n 1000000 > transcriptional.noise.rpm.txt

#R
Rscript /data/neurogen/pipeline/RNAseq/src/_fit.Tx.noise.R .
#for HCILB_SNDA: Dsig: 10**-1.105 == 0.079
Dsig=`tail -n1 transcriptional.noise.rpm.pvalues.txt`  

# step1: any regions with summit RPM > peakLevel and border > baseLevel
# =================
basalLevel=`tail -n1 $inputBG | cut -f2 -d'=' | cut -f1`
awk -vmin=$basalLevel '{OFS="\t"; if($4>=min) print $1,$2,$3,".",$4}' $inputBG | mergeBed -scores max > eRNA.tmp1

# step2: summit RPM >=Dsig (density with p<0.05)
# =================
awk -vD=$Dsig '{OFS="\t"; if($4>=D) print $1,$2,$3,".",$4}' eRNA.tmp1 | mergeBed -d 100 -scores max > eRNA.tmp2

# step3: located in non-generic regions (e.g. 500bp away from any annotated exons),
# =================
intersectBed -a eRNA.tmp2 -b ../toExclude.bed -v > eRNA.tmp3

# step4: length > 100nt
# =================
awk '{OFS="\t"; if(($3-$2)>100) print $1,$2,$3,$1"_"$2"_"$3}' eRNA.tmp3 > eRNA.tmp4

# step6: don't contain any splicing sites (donor or acceptor from trinity/cufflinks de novo assembly)
# =================
# cd ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/
# cat cufflinks-cuffmerge/merged.bed trinity-cuffmerge/all_strand_spliced.chr.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' | bed12ToBed6 | awk '{OFS="\t"; print $1, $2-10,$2+10; print $1,$3-10,$3+10;}' | sortBed | uniq > trinitycufflinks.merged.splicingsites.flanking20nt.bed

# more than 10 splicing reads in at least 5 samples
# for i in  ~/neurogen/rnaseq_PD/run_output/*/junctions.bed; do awk '{OFS="\t"; if($5>10) { split($11,a,","); split($12,b,","); print $1,$2+a[1]-10,$2+a[1]+10; print $1,$2+b[2]-10,$2+b[2]+10}}' $i | sortBed | uniq; done | sort | uniq -c | awk '{OFS="\t"; if($1>5) print $2,$3,$4}' > ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed
 
intersectBed -a eRNA.tmp4 -b ~/neurogen/rnaseq_PD/results2/merged/denovo_assembly/tophatjunctions.merged.splicingsites.flanking20nt.bed -v > eRNA.tmp5

# step5: calculate the significance of eRNA
# =================
#1: create 100,000 random regions (400bp each) as background and calculate their signals
for i in ~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw ~/neurogen/rnaseq_PD/results2/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bw;
do
    bsub -q normal -n 1 "bedtools shuffle -excl ../toExclude.bed -noOverlapping -i eRNA.tmp5 -g $GENOME/Annotation/Genes/ChromInfo.txt | bigWigAverageOverBed $i stdin stdout | cut -f1,5 > $i.rdbg";
    bsub -q normal -n 1 "bigWigAverageOverBed $i eRNA.tmp5 stdout | cut -f1,5 | sort -k1,1 > $i.eRNA.meanRPM"
done

### 2: distribution of random background, in order to define the cutoff with p=0.0001 significance
R
# significance
path=c("~/neurogen/rnaseq_PD/results2/merged/trimmedmean.uniq.normalized.HCILB_SNDA.bw", "~/neurogen/rnaseq_PD/run_output/[HI]*_SNDA*rep[0-9]/uniq/accepted_hits.normalized.bw")

## read in only the 90 HCILB subjects/samples (w/ genotype)
IDs=read.table('~/neurogen/rnaseq_PD/results2/merged/RNAseqID.wGenotyped.list',stringsAsFactors =F)[,1]
IDs=IDs[grep("^[HI].*_SNDA", IDs)]
EXP=data.frame(); PV=data.frame(); QV=data.frame(); id="locus"

pdf("background.RNAseq.cummulative.plot.pdf")
for(i in Sys.glob(path)){
    ii=ifelse(grepl("merged", i), sub(".*merged/(.*).bw.*","\\1", i), sub(".*run_output/(.*)/uniq.*","\\1", i));
    if(! (ii %in% IDs || grepl("trimmedmean",ii))) next;
    print(i)
    # read background
    df=read.table(paste(i,"rdbg",sep="."), header=F)[,2] # mean RPM (mean0 from bigWigAverageOverBed)
    Fn=ecdf(df)
    
    # plot the cummulative plot
    plot(Fn, verticals = TRUE, do.points = FALSE, main=ii, ylim=c(0.99, 1), xlab="average RPM", ylab="cummulative percentage (approx. 1-p)")
    inv_ecdf <- function(f){ x <- environment(f)$x; y <- environment(f)$y; approxfun(y, x)}; g <- inv_ecdf(Fn);
    abline(h=0.999, v=g(0.999), col='red', lty=2, lwd=1)
    points(g(0.999), 0.999, col='red', pch=19)
    text(g(0.999), 0.999, round(g(0.999),2), cex=5, adj=c(0,1))
    
    if(grepl("trimmedmean",ii)) next;
    id=c(id, ii)

    # read expression
    expression=read.table(paste(i,"eRNA.meanRPM",sep="."), header=F)
    pvalue=as.numeric(format(1-Fn(expression[,2]), digits=3));
    qvalue=as.numeric(format(p.adjust(pvalue, "BH"), digits=3));
    write.table(cbind(expression[,1:2], pvalue=pvalue, qvalue=qvalue), file=paste(i,"eRNA.meanRPM.significance",sep="."), quote=F, sep ="\t", col.names =F, row.names=F)
    
    # merge
    if(ncol(EXP)==0) { EXP=expression; expression[,2]=pvalue; PV=expression; expression[,2]=qvalue; QV=expression; }
    else {EXP=cbind(EXP, expression[,2]); PV=cbind(PV, pvalue); QV=cbind(QV, qvalue); }
}
dev.off()

colnames(EXP)=id; colnames(PV)=id; colnames(QV)=id;
rM=rowMeans(QV[,-1]<=0.05)
write.table(EXP[rM>0.25,], "eRNA.90samples.meanRPM.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(PV[rM>0.25,], "eRNA.90samples.pvalue.xls", col.names=T, row.names=F, sep="\t", quote=F)
write.table(QV[rM>0.25,], "eRNA.90samples.qvalue.xls", col.names=T, row.names=F, sep="\t", quote=F)

pdf("eRNA.90samples.qvalue.hist.pdf", width=8, height=6)
h=hist(rM, breaks=80, xlim=c(0,1), main="",xlab=expression("Percentage of HC/ILB SNDA samples (out of 90) with q-value" <= "0.05"), ylab="Count of HiTNEs", freq=T)
abline(v=0.250, lty=2, col='red')
legend('topright', c(bquote(.(sum(rM>0.25)) ~ "HiTNEs"), expression("with q-value" <= "0.05"), "in at least 25% of samples"),  bty='n', text.col='red', cex=1.5)
dev.off()

q('no')
## R end

awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) print a[1],a[2],a[3],$1}' eRNA.90samples.meanRPM.xls > eRNA.bed


## merge menaRPM for all samples
R
path="~/neurogen/rnaseq_PD/run_output/*/uniq/accepted_hits.normalized.bw"
EXP=data.frame(); id="locus"
for(i in Sys.glob(path)){
    ii=sub(".*run_output/(.*)/uniq.*","\\1", i);
    print(i)
    id=c(id, ii)
    expression=read.table(paste(i,"eRNA.meanRPM",sep="."), header=F)
    if(ncol(EXP)==0) {
      EXP=expression; 
    } else {EXP=cbind(EXP, expression[,2]);}
}
colnames(EXP)=id; 
write.table(EXP, "eRNA.140samples.meanRPM.xls", col.names=T, row.names=F, sep="\t", quote=F)