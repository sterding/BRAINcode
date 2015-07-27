## script to generate data to draw aggregation plot for HITNE

#============================================================
# get bin signal
#============================================================
cd ~/eRNAseq

# always take the middle position
awk '{OFS="\t"; mid=int(($2+$3)/2); print $1, mid-1000, mid+1000,$4;}' eRNA.bed > eRNA.1kbp.mid.1kbp.bed

# if overlapping with any DNase peak (from fetal brain, only 1/10 of regions), take the summit of the peak; otherwise, take the center point of HiTNEs. Then expand [-1k,+1k] from whatever the point determined above.
intersectBed -a eRNA.bed -b externalData/DNase/DNase.Roadmap.fBrain.narrowPeak -wao | sort -k4,4 -k15,15nr | awk '{if($4!=id) {id=$4;print;}}' | awk '{OFS="\t"; mid=($5==".")?int(($2+$3)/2):($6+$14); print $1, mid-1000, mid+1000,$4;}'> eRNA.1kbp.summit.or.mid.1kbp.bed

cat eRNA.bed | while read chr start end name score strand rest
do
    N=`expr $end - $start`;
    s=`bigWigSummary externalData/DNase/DNase.Roadmap.fBrain.pval.bigwig -udcDir=/tmp -type=mean $chr $start $end $N 2>/dev/null | sed 's/n\/a/0/g'`
    [[ $s == "" ]] && whichmax=`expr $N / 2`
    whichmax=`echo $s | awk 'BEGIN{max=-9999;}{for(i=1;i<=NF;i++) if($i>max) {which=i;max=$i}}END{print which}'`
    [[ $name == "" ]] && name="$chr_$start_$end";
    echo $chr `expr $start + $whichmax - 1000` `expr $start + $whichmax + 1000` $name;
done | sed 's/ /\t/g' > eRNA.1kbp.summit.or.mid.1kbp.bed

for i in externalData/*/*.bigwig;
do
    [ -e $i.eRNA.1kbp.summit.or.mid.1kbp.100bins ] || toBinRegionsOnBigwig.sh $i eRNA.1kbp.summit.or.mid.1kbp.bed 100 > $i.eRNA.1kbp.summit.or.mid.1kbp.100bins &
done

#============================================================
# draw aggregation plot (to run in R console)
#============================================================
Rscript ~/neurogen/pipeline/RNAseq/src/eRNA.aggPlot.R