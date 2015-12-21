####################################
# Pipeline for CAGE data Analysis
# Authors: Xianjun DOng
# Email: xdong@rics.bwh.harvard.edu
# Date: 11/21/2014
# Version: 0.0
####################################
#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: $HOME/neurogen/pipeline/RNAseq/CAGE.pipeline.sh /data/neurogen/CAGE_PDBrainMap/rawfiles 2"
  exit
fi


########################
## 0. setting 
########################
pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

# project folders
input_dir=$1  # input_dir=/data/neurogen/CAGE_PDBrainMap/rawfiles
batch=$2

# create the subfolders (e.g. processed, for_display, results)
processed=$input_dir/../processed
[ -d $processed ] || mkdir $processed

output_dir=$input_dir/../output_dir
[ -d $output_dir ] || mkdir $output_dir

fordisplay_dir=$input_dir/../for_display
[ -d $fordisplay_dir ] || mkdir $fordisplay_dir

result_dir=$input_dir/../results 
[ -d $result_dir ] || mkdir $result_dir

########################
# 0. sort fastq from barcode, clean up reads
########################

[ "$batch" == "1" ] && \
  cd $input_dir/batch1 && \ 
  ls separate_lane*.tar.gz | parallel tar -zxvf {} && \ 
  for i in lane6_Undetermined_L006_R1_0*.fastq; do echo $i; _CAGE.sorter.awk -v BARCODE='AGA|CTT|GAT' $i >> unsorted.lane6.fq; done && \ 
  for i in lane7_Undetermined_L007_R1_0*.fastq; do echo $i; _CAGE.sorter.awk -v BARCODE='ACA|ACT|ACG' $i >> unsorted.lane7.fa; done && \ 
  rm lane6_Undetermined_L006_R1_0*.fastq lane7_Undetermined_L007_R1_0*.fastq && \ 
  cat barcodes_lane* | awk '{split($1,a,"_");print "ln -fs sorted."$2".fq "a[1]".fastq";}' | bash;

[ "$batch" == "2" ] && \
  cd $input_dir/batch2 && \ 
  while read a b c rest; do echo ${a##*\/}; zcat ${a##*\/} | _CAGE.sorter.awk -v BARCODE="$c" >> unsorted.$b.fq; ln -fs sorted.$c.fq $b.fastq; done < barcodes.txt

cd $input_dir
ln -fs batch1/MCL7798.fastq ILB_MCL7798_SN_1_rep1.fastq
ln -fs batch1/MCL7878.fastq ILB_MCL7878_SN_1_rep1.fastq
ln -fs batch1/MD5247.fastq HC_MD5247_SN_1_rep1.fastq
ln -fs batch1/UW479.fastq HC_UWA479_SN_1_rep1.fastq # corrected
ln -fs batch1/UW616.fastq HC_UWA616_SN_1_rep1.fastq # corrected
ln -fs batch1/UW734.fastq PD_UWA734_SN_1_rep1.fastq # corrected

ln -fs batch2/UWA479.fastq HC_UWA479_SN_2_rep2.fastq
ln -fs batch2/UWA616.fastq HC_UWA616_SN_2_rep2.fastq
ln -fs batch2/UWA5247.fastq HC_MD5247_SN_2_rep2.fastq  # corrected

########################
## 1. Processing per sample 
########################
for R1 in *.fastq;
do
    samplename=${R1/.fastq*/}
    bsub -J $samplename -oo $output_dir/$samplename/_mapping.log -eo $output_dir/$samplename/_mapping.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] "tophat -o $output_dir/$samplename --keep-fasta-order -p 8 --read-mismatches 2 --segment-length 10 --max-multihits 10 --no-novel-juncs --no-coverage-search $GENOME/Sequence/Bowtie2Index/genome $R1"
done
# continue until the jobs are done

########################
# 2. convert bam to bigwig
########################

cd $output_dir
for i in `ls -d *_rep*`; do
    echo $i;
    [ -e $i/accepted_hits.bam ] && bsub -J _bam2bigwig -oo $output_dir/$i/_bam2bw.log -eo $output_dir/$i/_bam2bw.log -q $QUEUE -n $CPU -M $MEMORY -R rusage[mem=$MEMORY] _CAGE_bam2bigwig.sh $i/accepted_hits.bam;
done

# continue until the jobs are done

########################
# 3. computing trimmed mean of bedGraph
########################

unionBedGraphs -i `ls */*plus.normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); sum=0; for(i=a+1;i<=(c-a);i++) sum+=j[i];return sum/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > ../results/trimmedmean.plus.normalized.bg
bedGraphToBigWig ../results/trimmedmean.plus.normalized.bg $GENOME/Annotation/Genes/ChromInfo.txt ../results/trimmedmean.plus.normalized.bw

unionBedGraphs -i `ls */*minus.normalized.bedGraph` | awk 'function trimmedMean(v, p) {c=asort(v,j); a=int(c*p); sum=0; for(i=a+1;i<=(c-a);i++) sum+=j[i];return sum/(c-2*a);} {OFS="\t"; n=1; for(i=4;i<=NF;i++) S[n++]=$i; tm=trimmedMean(S, 0.1); if(tm!=0) print $1,$2,$3,tm; s+=tm*($3-$2);}END{TOTAL=3095677412; bc=s/TOTAL; print "#basalCoverage="bc, "#"s"/"TOTAL;}' | awk '{OFS="\t"; if(id!=$4 || e!=$2 || chr!=$1) {if(chr!="") print chr,s,e,id; chr=$1;s=$2;e=$3;id=$4;} else {e=$3;}}END{print chr,s,e,id}' > ../results/trimmedmean.minus.normalized.bg
bedGraphToBigWig ../results/trimmedmean.minus.normalized.bg $GENOME/Annotation/Genes/ChromInfo.txt ../results/trimmedmean.minus.normalized.bw

for i in */*bw; do scp $i xd010@panda.dipr.partners.org:~/public_html/cage/version3/${i/\//.}; done
for i in ../results/*bw; do scp $i xd010@panda.dipr.partners.org:~/public_html/cage/version3/; done

########################
# 4. call bidirectional loci accoroding to Andersson et al 2014 Nature (https://github.com/anderssonrobin/enhancers)
########################
# prepare bed file for transcription initiation sites (5' bed6 files with score column quantifying the number of mapped reads)
for i in */accepted_hits.bed; do echo $i; sort -k1,1 -k2,2n -k6,6 $i | groupBy -g 1,2,3,6 -c 6 -o count | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2,$5,$4;}' > $i.merged.BED & done

awk '{print $0 >> $1}' accepted_hits.bed
for i in chr*; do echo $i; sort -k1,1 -k2,2n -k6,6 $i | groupBy -g 1,2,3,6 -c 6 -o count | awk '{OFS="\t"; print $1,$2,$3,$1"_"$2,$5,$4;}' > $i.merged.BED & done
cat chr*.merged.BED > 
ls -1 */accepted_hits.bed.merged.BED > ctss.filelist
~/bin/enhancers/scripts/bidir_enhancers -s braincode -m ~/bin/enhancers/mask/hg19_neg_filter_500_merged.bed -f ctss.filelist -o .



