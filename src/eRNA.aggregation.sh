## script to generate data to draw aggregation plot for HITNE
## Note: update on Jun 7th, 2017 (see old archieve in github if needed)

#============================================================
# extract DNase peaks from Roadmap brain samples
#============================================================

cd ~/eRNAseq/externalData/DNase

## merge the uniform signal of individual samples into the final "DNase bigwig track"
grep DNase macs2signal.list | cut -f2 | awk '{print "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"$1}' | xargs -n 1 -P 8 wget -b
mkdir download; mv E* download/
bsub -q big -n 1 -M 10000 bigWigMerge download/E*DNase.pval.signal.bigwig merged.DNase.pval.signal.bg
bsub -q big -n 1 -M 10000 bedGraphToBigWig merged.DNase.pval.signal.bg $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size merged.DNase.pval.signal.bigwig
# scp to http://panda.partners.org/~xd010/tracks/

## define DNaes peaks 
## using the 638304 enhancer-being DNase peaks defined in 10 Roadmap "brain" samples (http://egg2.wustl.edu/roadmap/web_portal/DNase_reg.html#delieation) [SELECTED as final solution]
# download regions_enh_*.bed.gz for 10 brain samples (E067-074, E081 and E082)
zcat regions_enh_*.bed.gz | sortBed | mergeBed -i - | awk '{OFS="\t"; print $0,"peak"NR }'> regions_enh_merged.brain.bed
split -l 1000 regions_enh_merged.brain.bed tmp.regions_enh_merged.brain.bed
for i in tmp.regions_enh_merged.brain.bed*; do echo $i; bsub -q short -n 1 -M 1000 bash ~/pipeline/bin/bed2narrowpeak.sh $i merged.DNase.pval.signal.bigwig; done
cat tmp.regions_enh_merged.brain.bed*narrowPeak > regions_enh_merged.brain.narrowPeak
rm tmp.regions_enh_merged.brain.bed*

#============================================================
# get bin signal
#============================================================
cd ~/eRNAseq/HCILB_SNDA

## if overlapping with any DNase peak (from fetal brain, only 1/10 of regions), take the summit of the peak; otherwise, take the center point of HiTNEs. Then expand [-1k,+1k] from whatever the point determined above.
#intersectBed -a eRNA.bed -b ../externalData/DNase/DNase.Roadmap.fBrain.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed
# change to use the narrowPeak based on merged Roadmap brain 
intersectBed -a eRNA.bed -b ../externalData/DNase/regions_enh_merged.brain.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed
intersectBed -a eRNA.bed -b ../externalData/DNase/regions_enh_merged.brain.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-2000, mid+2000,$4;}'> eRNA.2kbp.summit.or.mid.2kbp.bed

## below code is similar as bin/bed2narrowPeak.sh + above line
# cat eRNA.bed | while read chr start end name score strand rest
# do
#     N=`expr $end - $start`;
#     s=`bigWigSummary externalData/DNase/DNase.Roadmap.fBrain.pval.bigwig -udcDir=/tmp -type=mean $chr $start $end $N 2>/dev/null | sed 's/n\/a/0/g'`
#     [[ $s == "" ]] && whichmax=`expr $N / 2`
#     whichmax=`echo $s | awk 'BEGIN{max=-9999;}{for(i=1;i<=NF;i++) if($i>max) {which=i;max=$i}}END{print which}'`
#     [[ $name == "" ]] && name="$chr_$start_$end";
#     echo $chr `expr $start + $whichmax - 1000` `expr $start + $whichmax + 1000` $name;
# done | sed 's/ /\t/g' > eRNA.1kbp.summit.or.mid.1kbp.bed

for i in ../externalData/*/*.bigwig;
do
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || bsub -q short -n 1 -M 1000 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins"
done

# for externalData/Histone/Roadmap -- too large to run bigWigSummmary as a whole
for i in ../externalData/Histone/Roadmap/*chr*bigwig;
do
    j=${i/*Merged./}
    chr=${j/.bigwig/}
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins.$chr ] || bsub -q short -n 1 -M 500 "toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed.$chr 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins.$chr"
done
for i in H3K4me1 H3K4me3 H3K27ac; do cat ../externalData/Histone/Roadmap/RoadmapBrain${i}Merged.chr*eRNA.1kbp.summit.or.mid.1kbp.100bins.chr* > ../externalData/Histone/$i.Roadmap.brainMerged.bigwig.eRNA.1kbp.summit.or.mid.1kbp.100bins; done
rm ../externalData/Histone/Roadmap/RoadmapBrain*Merged.chr*eRNA.1kbp.summit.or.mid.1kbp.100bins.chr*

awk '{print $0 >> "eRNA.2kbp.summit.or.mid.2kbp.bed."$1}' eRNA.2kbp.summit.or.mid.2kbp.bed
for i in ../externalData/Histone/Roadmap/*chr*bigwig;
do
    j=${i/*Merged./}
    chr=${j/.bigwig/}
    [ -e $i.eRNA.2kbp.summit.or.mid.2kbp.100bins.$chr ] || bsub -q short -n 1 -M 500 "toBinRegionsOnBigwig.sh $i eRNA.2kbp.summit.or.mid.2kbp.bed.$chr 100 > $i.eRNA.2kbp.summit.or.mid.2kbp.100bins.$chr"
done
for i in H3K4me1 H3K4me3 H3K27ac; do cat ../externalData/Histone/Roadmap/RoadmapBrain${i}Merged.chr*eRNA.2kbp.summit.or.mid.2kbp.100bins.chr* > ../externalData/Histone/$i.Roadmap.brainMerged2k.bigwig.eRNA.1kbp.summit.or.mid.1kbp.100bins; done
rm ../externalData/Histone/Roadmap/RoadmapBrain*Merged.chr*eRNA.2kbp.summit.or.mid.2kbp.100bins.chr*



#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.aggPlot.R