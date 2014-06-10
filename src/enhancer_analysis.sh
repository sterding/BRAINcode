## script to analysis enhancer regions

ANNOTATION=$GENOME/Annotation/Genes

cd ~projects/PD/results/eRAN/externalData

# regions to exclude (e.g. 500nt flanking of TSS, 200nt flanking of exons)
cat <(cut -f1-6 $ANNOTATION/gencode.v19.annotation.bed12 | flankBed -l 500 -r 0 -s -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai) <(bedtools bed12tobed6 -i $ANNOTATION/gencode.v19.annotation.bed12 | slopBed -g $GENOME/Sequence/WholeGenomeFasta/genome.fa.fai -b 200) > toExclude.exontss.bed

#============================================================
# Enhancers defined by ENCODE histone modifications
# e.g. active enhancers = H3K4me1 + H3k27ac - H3K27me3
# silent enhancers = H3K4me1 - H3k27ac + H3K27me3
#============================================================
cd ENCODE
url=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone
for i in K562 Gm12878 H1hesc;
do
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz;
    ##wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k04me1StdPk.broadPeak.gz;
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz;
    #wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz;
     
    intersectBed -a <(zcat wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz) -b <(zcat wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz) -wa -u | intersectBed -a - -b <(zcat wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz) -v -wa | awk -vi=$i '{OFS="\t"; $4=i; print}' > $i.silent.enhancer.broadPeak
    intersectBed -a <(zcat wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz) -b <(zcat wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz) -wa -u | intersectBed -a - -b <(zcat wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz) -v -wa | awk -vi=$i '{OFS="\t"; $4=i; print}' > $i.active.enhancer.broadPeak
    echo $i;
done

# P300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsSydhK562Gata1UcdUniPk.narrowPeak.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsBroadK562P300UniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "K562";}' | intersectBed -a - -b K562.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.exontss.bed -v -wa > K562.active.enhancer.broadPeak.distal.wP300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibGm12878P300Pcr1xUniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "Gm12878";}' | intersectBed -a - -b Gm12878.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.exontss.bed -v -wa > Gm12878.active.enhancer.broadPeak.distal.wP300
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibH1hescP300V0416102UniPk.narrowPeak.gz | zcat |awk '{OFS="\t"; peak=$2+$10; print $1, peak-500, peak+500, "H1hesc";}' | intersectBed -a - -b H1hesc.active.enhancer.broadPeak -u | intersectBed -a - -b ../toExclude.exontss.bed -v -wa > H1hesc.active.enhancer.broadPeak.distal.wP300

# merge
sort -k1,1 -k2,2n *active.enhancer.broadPeak.distal.wP300 > ENCODE.distal.bed

# bigwig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k04me1StdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k27acStdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k27me3StdSig.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneHelas3H3k4me3StdSig.bigWig

ln -s wgEncodeBroadHistoneHelas3H3k04me1StdSig.bigWig ENCODE.H3k4me1.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k27acStdSig.bigWig ENCODE.H3k27ac.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k27me3StdSig.bigWig ENCODE.H3k27me3.bigwig
ln -s wgEncodeBroadHistoneHelas3H3k4me3StdSig.bigWig ENCODE.H3k4me3.bigwig

#============================================================
# Enhancers defined by TFBS HOT region
#============================================================
cd ../HOT
curl -s http://stanford.edu/~claraya/metrn/data/hot/regions/hs/maphot_hs_selection_reg_cx_occP05_any.bed | cut -f1-6  > maphot_hs_selection_reg_cx_occP05_any.bed6
# distal HOT regions 
intersectBed -a maphot_hs_selection_reg_cx_occP05_any.bed6 -b ../toExclude.exontss.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > HOT.distal.bed

# clustered TFBS (count all Peaks for 161 transcription factors in 91 cell types)
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
gunzip -c wgEncodeRegTfbsClusteredV3.bed.gz | awk '{OFS="\t"; print $1,$2,$3,$4,0,"+",$2,$3,"0",$6,$7,$8}' > wgEncodeRegTfbsClusteredV3.bed12
sort -k1,1 wgEncodeRegTfbsClusteredV3.bed12 | bedItemOverlapCount hg19 -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | sed 's/ /\t/g' > wgEncodeRegTfbsClusteredV3.bg
bedGraphToBigWig wgEncodeRegTfbsClusteredV3.bg $ANNOTATION/ChromInfo.txt wgEncodeRegTfbsClusteredV3.bw
ln -fs wgEncodeRegTfbsClusteredV3.bw HOT.bigwig

#============================================================
# Enhancers defined by FANTOM5 CAGE
#============================================================
cd ../CAGE
curl -s http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed | awk '{OFS="\t"; peak=$8; if(peak>500) print $1, peak-500, peak+500;}' > CAGE.distal.bed

wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.cell_line.hCAGE/epitheloid%20carcinoma%20cell%20line:%20HelaS3%20ENCODE,%20biol_rep1.CNhs12325.10815-111B5.hg19.ctss.fwd.bw
wget http://fantom.gsc.riken.jp/5/datahub/hg19/ctss/human.cell_line.hCAGE/epitheloid%20carcinoma%20cell%20line:%20HelaS3%20ENCODE,%20biol_rep1.CNhs12325.10815-111B5.hg19.ctss.rev.bw

ln -fs epitheloid*ctss.fwd.bw CAGE.fwd.bigwig
ln -fs epitheloid*ctss.rev.bw CAGE.rev.bigwig

#============================================================
# Enhancers defined by conservation
#============================================================
cd ../Conservation
# Ancora HCNEs
curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_hg19_danRer7_70pc_50col.bed.gz | gunzip | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp
ln -sf HCNE_hg19_danRer7_70pc_50col.bed.trimedto1000bp Conservation.distal.bed


## bigwig file
#curl -s http://ancora.genereg.net/downloads/hg19/vs_zebrafish/HCNE_density_hg19_danRer7_70pc_50col.wig.gz | gunzip | wigToBigWig stdin $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size HCNE_density_hg19_danRer7_70pc_50col.bw
#ln -s HCNE_density_hg19_danRer7_70pc_50col.bw HCNE.distal.bigwig

# phyloP
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate ./
zcat vertebrate/chr*.gz | gzip -c > vertebrate/phyloP46way.wigFix.gz
wigToBigWig vertebrate/phyloP46way.wigFix.gz $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size vertebrate/phyloP46way.wigFix.bigwig

ln -sf vertebrate/phyloP46way.wigFix.bigwig Conservation.bigwig

## 29 mammals conservation (in hg18)
#curl -s http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz | gunzip > 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt
#intersectBed -a way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt -b ../toExclude.exontss.bed -v -wa
#
## bigwig file
#wget http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/omega.12mers.wig.gz
#wigToBigWig omega.12mers.wig $GENOME/Sequence/WholeGenomeFasta/hg18.chrom.size omega.12mers.hg18.bigwig

#============================================================
# Enhancers defined by RNAseq
#============================================================
# uniq mapping RNAseq with (1) RPM>0.5 (2) non-genic and (3) length>=50   --> uniq/accepted_hits.normalized.eRNA.bed in each sample folder
cd ../RNAseq

exons=$ANNOTATION/gencode.v13.annotation.gtf.exons.bed
for i in /data/neurogen/rnaseq_PD/results/merged/*trimmedmean.uniq.normalized.bedGraph; do echo $i; awk '$4>0.05' $i | mergeBed | intersectBed -a - -b $exons -v | awk '($3-$2)>=50' > ${i/bedGraph/eRNA.bed} & done
intersectBed -a /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.eRNA.bed -b ../toExclude.exontss.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > RNAseq.distal.bed

ln -s /data/neurogen/rnaseq_PD/results/merged/HC_SNDA.trimmedmean.uniq.normalized.bw RNAseq.distal.bigwig
#============================================================
# Enhancers defined by DNase
#============================================================
cd ../DNase
# ref: http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=377332155_oXR9pqyaZLzzFHxgO4t0YzyQseGN&c=chr1&g=wgEncodeRegDnaseClusteredV2
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV2.bed.gz
zcat wgEncodeRegDnaseClusteredV2.bed.gz | awk '$4>20 && $5>800' | intersectBed -a - -b ../toExclude.exontss.bed -v | awk '{OFS="\t"; mid=int(($3+$2)/2); if(mid>500) print $1, mid-500, mid+500;}' > wgEncodeRegDnaseClusteredV2.bed  # at least occuring in 20 cell lines (out of 125), and 800 cluster score (out of 1000) --> 100180 regions remained

ln -s wgEncodeRegDnaseClusteredV2.bed DNase.distal.bed


#wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz
#wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseHelas3PkRep1.narrowPeak.gz
#wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseGm12878PkRep1.narrowPeak.gz
#wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseH1hescPkRep1.narrowPeak.gz
#zcat wgEncodeUwDnaseH1hescPkRep1.narrowPeak.gz | intersectBed -a - -b ../toExclude.exontss.bed -v | cut -f1-3 > DNase.Helas3.bed


wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseHelas3RawRep1.bigWig
wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562RawRep1.bigWig
wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseH1hescRawRep1.bigWig
wb http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseGm12878RawRep1.bigWig

ln -s wgEncodeUwDnaseHelas3RawRep1.bigWig DNase.Helas3.bigwig

#============================================================
# get bin signal
#============================================================
cd ~/neurogen/rnaseq_PD/results/eRNA/externalData

for i in ENCODE CAGE HOT Conservation DNase RNAseq;
do
    for j in */*bigwig; do
        echo $i, $j;
        [ -f $j.$i.bins ] || toBinRegionsOnBigwig.sh $j $i/$i.distal.bed 200 > $j.$i.bins &
    done
    #for expressed regions only (for DoD grant)
    intersectBed -a $i/$i.distal.bed -b RNAseq/RNAseq.distal.bed -u > $i/${i}_n_RNAseq.distal.bed
    [ -f RNAseq/RNAseq.${i}_n_RNAseq.bins ] || toBinRegionsOnBigwig.sh RNAseq/RNAseq.bigwig $i/${i}_n_RNAseq.distal.bed 200 > RNAseq/RNAseq.${i}_n_RNAseq.bins &
done


#============================================================
# draw aggregation plot (to run in R console)
#============================================================
output_dir="~/neurogen/rnaseq_PD/results/eRNA/externalData";
pdf("aggregation.enhancers.pdf", width=18, height=16)
par(mfcol=c(7,6))
for(mark in c("CAGE","ENCODE","HOT","Conservation","DNase","RNAseq")){
    print(mark);
    
    n=nrow(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F))
    
    #CAGE
    df1=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    df2=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.rev.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    o=par(mar = c(2,4,4,2))
    plot(df1, type="l", col="green", ylab="Mean signal of CAGE", xlab="", xaxt="n", main=paste(mark, "-defined enhancers\n(N=",formatC(n,format="d",big.mark=","),")",sep=""), ylim=range(df1,df2)) #c(-max(df1,df2),max(df1,df2)))
    points(df2, type="l", col="blue")
    legend("topright", c("CAGE +","CAGE -"), col=c("green","blue"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #ENCODE
    df1=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27ac.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df2=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me1.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df3=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df4=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    par(mar = c(2,4,2,2))
    plot(df1, type="l", col="green", ylab="Mean signal of histone marks", xlab="", xaxt="n", ylim=range(df1,df2, df3, df4, 10)) 
    points(df2, type="l", col="blue")
    points(df3, type="l", col="gold")
    points(df4, type="l", col="black")
    legend("topright", c("H3K27ac","H3K4me1", "H3K4me3","H3K27me3"), col=c("green","blue", "gold","black"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #HOT
    df=apply(read.table(paste(output_dir, "HOT", paste("HOT.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean TFBS occurances", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #Conservation 
    df=apply(read.table(paste(output_dir, "Conservation", paste("Conservation.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean phyloP score", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #DNase
    df=apply(read.table(paste(output_dir, "DNase", paste("DNase.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean DNase peak occurances", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #RNAseq
    df=apply(read.table(paste(output_dir, "RNAseq", paste("RNAseq.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    par(mar = c(4,4,2,2))
    plot(df, type="l", col="green", ylab="Mean signal of RNAseq", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))

    # intersect with RNAseq
    df=read.table(paste(output_dir, "RNAseq", paste("RNAseq",paste(mark,"_n_RNAseq",sep=""),"bins",sep="."), sep="/"), header=F)
    n=nrow(df)
    df=apply(df[,-1], 2, function(x) mean(x, trim=0.01))
    par(mar = c(4,4,4,2))
    plot(df, type="l", col="red", ylab="Mean signal of RNAseq", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n", main=paste(mark, " enhancers expressed in SNDA\n(N=",formatC(n,format="d",big.mark=","),")",sep="")) 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
}
dev.off()


#============================================================
# draw plot for DoD grant  (run in R console)
#============================================================
output_dir="~/neurogen/rnaseq_PD/results/eRNA/externalData";
pdf("aggregation.enhancers.forDoD.pdf", width=6, height=12)
par(mfcol=c(6,2), mar = c(5,4,4,2))
for(mark in c("CAGE","ENCODE")){
    print(mark);
    
    n=nrow(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F))
    
    #CAGE
    df1=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.fwd.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    df2=apply(read.table(paste(output_dir, "CAGE", paste("CAGE.rev.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=10/n))
    plot(df1, type="l", col="green", ylab="Franction of CAGE tag", xlab="", xaxt="n", main=paste(mark, "-defined enhancers\n(N=",formatC(n,format="d",big.mark=","),")",sep=""), ylim=range(df1,df2)) #c(-max(df1,df2),max(df1,df2)))
    points(df2, type="l", col="blue")
    legend("topright", c("CAGE +","CAGE -"), col=c("green","blue"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #ENCODE
    df1=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27ac.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df2=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me1.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df3=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k4me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    df4=apply(read.table(paste(output_dir, "ENCODE", paste("ENCODE.H3k27me3.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df1, type="l", col="green", ylab="Mean signal of histone marks", xlab="", xaxt="n", ylim=range(df1,df2, df3, df4)) 
    points(df2, type="l", col="blue")
    points(df3, type="l", col="gold")
    points(df4, type="l", col="black")
    legend("topright", c("H3K27ac","H3K4me1", "H3K4me3","H3K27me3"), col=c("green","blue", "gold","black"), lty=1, bty='n')
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #HOT
    df=apply(read.table(paste(output_dir, "HOT", paste("HOT.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean TFBS occurances", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #Conservation 
    df=apply(read.table(paste(output_dir, "Conservation", paste("Conservation.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean phyloP score", xlab="", xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #DNase
    df=apply(read.table(paste(output_dir, "DNase", paste("DNase.bigwig",mark,"bins",sep="."), sep="/"), header=F)[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="green", ylab="Mean DNase peak occurances", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n") 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
    
    #RNAseq
    df=read.table(paste(output_dir, "RNAseq", paste("RNAseq.bigwig",paste(mark,"_n_RNAseq",sep=""),"bins",sep="."), sep="/"), header=F)
    n=nrow(df)
    df=apply(df[,-1], 2, function(x) mean(x, trim=0.01))
    plot(df, type="l", col="red", ylab="Mean RNAseq signal", xlab=paste("Center position of\n",mark, "-defined enhancers",sep=""), xaxt="n", main=paste(mark, " enhancers expressed in SNDA\n(N=",formatC(n,format="d",big.mark=","),")",sep="")) 
    axis(1, at=50*c(0:4), labels=c("-500","","0","","500"))
}
        
dev.off()
