# Uage:
# bash ~/neurogen/pipeline/RNAseq/modules/_CAGE_make_trackDb.sh > ~/neurogen/CAGE_PDBrainMap/for_display/trackDb.CAGE.txt
# scp ~/neurogen/CAGE_PDBrainMap/for_display/trackDb.CAGE.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
# to get the sample list: ls -1 ~/neurogen/cage_PD/run_output/ | sed 's/_/\t/g' | cut -f2 | sort -u | awk '{printf $1"="$1" "}'
#
# color selected from: http://colorbrewer2.org/

# output the general header [superTrack]

echo "track CAGE_BRAINCODE
shortLabel CAGE
longLabel BRAINCODE CAGE (version 1)
dataVersion Version 1 (Nov 2014)
type bed 3
visibility full
boxedCfg on
priority 24
superTrack on show
subGroup1 signal Signal RPM=RPM raw=Raw

"
# output the individual track [container]

for i in ~/neurogen/CAGE_PDBrainMap/processed/*.plus.bw; do
    i=${i/*\//}
    samplename=${i/.plus*/}
    
echo "
    track ${samplename}_multiwig
    container multiWig
    configurable on
    shortLabel raw $samplename
    longLabel CAGE tags count for $samplename
    aggregate transparentOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 64:64:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE
    subGroups signal=raw
    
        track ${samplename}_multiwig_fwd
        shortLabel CAGE tags count for $samplename (fwd)
        longLabel CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/$samplename.plus.bw
        color 255,0,0
        parent ${samplename}_multiwig
        
        
        track ${samplename}_multiwig_rev
        shortLabel CAGE tags count for $samplename (rev)
        longLabel CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/$samplename.minus.bw
        color 0,0,255
        parent ${samplename}_multiwig
        
"
done


# output the individual track [container] -- normalized

for i in ~/neurogen/CAGE_PDBrainMap/processed/*.plus.normalized.bw; do
    i=${i/*\//}
    samplename=${i/.plus*/}
    
echo "
    track ${samplename}_multiwignormalized
    container multiWig
    configurable on
    shortLabel rpm $samplename
    longLabel normalized CAGE tags count for $samplename
    aggregate transparentOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 64:64:11
    showSubtrackColorOnUi on
    priority 1.4
    parent CAGE_BRAINCODE
    subGroups signal=RPM
    
        track ${samplename}_multiwignormalized_fwd
        shortLabel normalized CAGE tags count for $samplename (fwd)
        longLabel normalized CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/$samplename.plus.normalized.bw
        color 255,0,0
        parent ${samplename}_multiwignormalized
        
        
        track ${samplename}_multiwignormalized_rev
        shortLabel normalized CAGE tags count for $samplename (rev)
        longLabel normalized CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/$samplename.minus.normalized.bw
        color 0,0,255
        parent ${samplename}_multiwignormalized
"
done