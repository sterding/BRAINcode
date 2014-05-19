## script to analysis enhancer regions

cd /Users/xdong/projects/PD/results/eRAN/externalData

#============================================================
# Enhancers defined by ENCODE histone modifications
# e.g. active enhancers = H3K4me1 + H3k27ac - H3K27me3
# silent enhancers = H3K4me1 - H3k27ac + H3K27me3
#============================================================
cd ENCODE
for i in K562 Gm12878 H1hesc;
do
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz;
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz;
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz;
    
    intersectBed -a <(gzcat wgEncodeBroadHistone${i}H3k4me1StdPk.broadPeak.gz) -b <(gzcat wgEncodeBroadHistone${i}H3k27me3StdPk.broadPeak.gz) -wa -u | intersectBed -a - -b <(gzcat wgEncodeBroadHistone${i}H3k27acStdPk.broadPeak.gz) -v -wa > $i.silent.enhancer.broadPeak
    

#============================================================
# Enhancers defined by TFBS HOT region
#============================================================

#============================================================
# Enhancers defined by FANTOM5 CAGE
#============================================================


#============================================================
# Enhancers defined by conservation
#============================================================

#============================================================
# Enhancers defined by 
#============================================================

