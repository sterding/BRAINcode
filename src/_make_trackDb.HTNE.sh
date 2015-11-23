# Uage:
# cd ~/neurogen/rnaseq_PD/for_display/; _make_trackDb.HTNE.sh ; cd -
#
# color selected from: http://colorbrewer2.org/
# Note: use the Google doc to define/change the color code: https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/edit#gid=1995457670

curl -sk "https://docs.google.com/spreadsheets/d/1Sp_QLRjFPW6NhrjNDKu213keD_H9eCkE16o7Y1m35Rs/pub?gid=1995457670&single=true&output=tsv" > /tmp/braincode.colors.txt

echo "# Making merged bigBed files for major and minor groups"
###############################################################

for i in HCILB_SNDA HC_nonNeuron HC_PY HC_FB HC_PBMC ILB_SNDA PD_SNDA HC_MCPY HC_TCPY; 
do 
  RGB=`fgrep -w $i /tmp/braincode.colors.txt | cut -f3`
  awk -vRGB=$RGB -vI=$i '{OFS="\t";print $1,$2,$3,$4"."I,1000,".",$2,$3,RGB}' ~/eRNAseq/$i/eRNA.fdr.bed | sortBed > HTNE.fdr.$i.bed
  bedToBigBed HTNE.fdr.$i.bed $GENOME/Annotation/Genes/ChromInfo.txt HTNE.fdr.$i.bb -tab -type=bed9
done

chmod 644 HTNE*.bb; scp HTNE*.bb xd010@panda.dipr.partners.org:~/public_html/rnaseq_PD/version3/merged

echo "# Making trackDb for UCSC"
###############################################################

echo "track HTNE
shortLabel HTNE
longLabel BRAINCODE HTNE (corrected with fdr <0.05)
dataVersion Version 4 (Nov 2015)
type bed 3
visibility full
boxedCfg on
priority 3
superTrack on show

    track myComposite
    compositeTrack on
    parent HTNE
    shortLabel HTNE
    longLabel BRAINCODE HTNE (corrected with fdr <0.05)
    visibility dense
    allButtonPair on
" > trackDb.HTNE.txt;

CELLS=(HC_FB HC_PBMC HC_MCPY HC_TCPY PD_SNDA ILB_SNDA HCILB_SNDA HC_PY HC_nonNeuron)
tLen=${#CELLS[@]}

for I in `seq 1 $tLen`;
do
  i=${CELLS[`expr $I - 1`]}  
  echo "
        track HTNE_$i
        shortLabel $i
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version3/merged/HTNE.fdr.$i.bb
        priority 3.$I
        type bigBed 9 .
        itemRgb on
        visibility dense
        parent myComposite on
  " >> trackDb.HTNE.txt;
done

chmod 644 trackDb.HTNE.txt; scp trackDb.HTNE.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/

echo "# DONE"
