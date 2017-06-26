# Uage:
# bash ~/neurogen/pipeline/RNAseq/modules/_CAGE_make_trackDb.sh > ~/neurogen/CAGE_PDBrainMap/for_display/trackDb.CAGE.v3.txt
# chmod 644 ~/neurogen/CAGE_PDBrainMap/for_display/trackDb.CAGE.v3.txt
# rsync -azv ~/neurogen/CAGE_PDBrainMap/for_display/trackDb.CAGE.v3.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
# to get the sample list: ls -1 ~/neurogen/cage_PD/run_output/ | sed 's/_/\t/g' | cut -f2 | sort -u | awk '{printf $1"="$1" "}'
#
# color selected from: http://colorbrewer2.org/

# output the general header [superTrack]

echo "track CAGE_BRAINCODE
shortLabel CAGE
longLabel BRAINCODE CAGE (version 3)
dataVersion Version 3 (Dec 2015)
type bed 3
visibility full
boxedCfg on
priority 24
superTrack on show
"

# output the combined track [container multiWig]

echo "track merged_by_trimmedmean_CAGE
    shortLabel CAGE_combined
    longLabel BRAINCODE CAGE combined tracks
    configurable on
    container multiWig
    aggregate solidOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 100:64:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE
    
        track merged_multiwig_fwd
        shortLabel CAGE tags count for merged samplename (fwd)
        longLabel CAGE tags count for merged samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/trimmedmean.plus.normalized.bw
        color 255,0,0
        parent merged_by_trimmedmean_CAGE
        
        
        track merged_multiwig_rev
        shortLabel CAGE tags count for merged samplename (rev)
        longLabel CAGE tags count for merged samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/trimmedmean.minus.normalized.bw
        color 0,0,255
        altColor 0,0,255
        parent merged_by_trimmedmean_CAGE
"

echo "track	HCILB_merged_CAGE
    shortLabel Total count (BC, SN)
    longLabel BRAINCODE HCILB_merged CAGE tracks
    configurable on
    container multiWig
    aggregate solidOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 100:64:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE
    
        track HCILB_merged_CAGE_fwd
        shortLabel CAGE tags count for merged samplename (fwd)
        longLabel CAGE tags count for merged samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/HCILB_merged_accepted_hits.plus.bw
        color 255,0,0
        parent HCILB_merged_CAGE
        
        
        track HCILB_merged_CAGE_rev
        shortLabel CAGE tags count for merged HCILB samples (rev)
        longLabel CAGE tags count for merged HCILB samples reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/HCILB_merged_accepted_hits.minus.bw
        color 0,0,255
        altColor 0,0,255
        parent HCILB_merged_CAGE
"

# add FANTOM5 total count
echo "track	FANTOM_total_count
    shortLabel Total count (FANTOM)
    longLabel FANTOM Total counts of CAGE tracks
    configurable on
    container multiWig
    aggregate solidOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 100:64:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE
    
          track FANTOM_total_count_Fwd
          shortLabel Total counts of CAGE reads (fwd)
          longLabel Total counts of CAGE reads forward
          type bigWig
          negateValues off
          description Total Pooled Counts(Fwd)
          bigDataUrl http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.fwd.bw
          color 255,0,0
          altColor 255,0,0
          parent FANTOM_total_count

          track FANTOM_total_count_Rev
          shortLabel Total counts of CAGE reads (rev)
          longLabel Total counts of CAGE reads reverse
          type bigWig
          negateValues on
          description Total Pooled Counts(Rev)
          bigDataUrl http://fantom.gsc.riken.jp/5/datahub/hg19/reads/ctssTotalCounts.rev.bw
          color 0,0,255
          altColor 0,0,255
          parent FANTOM_total_count
"


# output the individual track [container]

for i in ~/neurogen/CAGE_PDBrainMap/output_dir/*/*.plus.bw; do
    i=${i/*output_dir\//}
    samplename=${i/\/*/}
    
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
    visibility hide
    maxHeightPixels 64:50:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE
    subGroups signal=raw
    
        track ${samplename}_multiwig_fwd
        shortLabel CAGE tags count for $samplename (fwd)
        longLabel CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.accepted_hits.plus.bw
        color 255,0,0
        parent ${samplename}_multiwig
        
        
        track ${samplename}_multiwig_rev
        shortLabel CAGE tags count for $samplename (rev)
        longLabel CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.accepted_hits.minus.bw
        color 0,0,255
        altColor 0,0,255
        parent ${samplename}_multiwig
        
"
done

# output the individual track [container] -- normalized

for i in ~/neurogen/CAGE_PDBrainMap/output_dir/*/*plus.normalized.bw; do
    i=${i/*output_dir\//}
    samplename=${i/\/*/}
    
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
    maxHeightPixels 64:50:11
    showSubtrackColorOnUi on
    priority 1.4
    parent CAGE_BRAINCODE
    subGroups signal=RPM
    
        track ${samplename}_multiwignormalized_fwd
        shortLabel normalized CAGE tags count for $samplename (fwd)
        longLabel normalized CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.accepted_hits.plus.normalized.bw
        color 255,0,0
        parent ${samplename}_multiwignormalized
        
        
        track ${samplename}_multiwignormalized_rev
        shortLabel normalized CAGE tags count for $samplename (rev)
        longLabel normalized CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.accepted_hits.minus.normalized.bw
        color 0,0,255
        altColor 0,0,255
        parent ${samplename}_multiwignormalized
"
done

############################################################
## -------------------- DZNE alignment ------------------ ##
############################################################

for i in ~/neurogen/CAGE_PDBrainMap/data_upload_scherzer_20151124/bam_files/*.plus.bw; do
    i=${i/*\//}
    samplename=${i/.plus*/}
    
echo "
    track ${samplename}_multiwig_DZNE
    container multiWig
    configurable on
    shortLabel DZNE raw $samplename
    longLabel DZNE CAGE tags count for $samplename
    aggregate transparentOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility hide
    maxHeightPixels 64:50:11
    showSubtrackColorOnUi on
    priority 1.3
    parent CAGE_BRAINCODE

        track ${samplename}_multiwig_DZNE_fwd
        shortLabel DZNE CAGE tags count for $samplename (fwd)
        longLabel DZNE CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.plus.bw
        color 255,0,0
        parent ${samplename}_multiwig_DZNE
        
        
        track ${samplename}_multiwig_DZNE_rev
        shortLabel DZNE CAGE tags count for $samplename (rev)
        longLabel DZNE CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.minus.bw
        color 0,0,255
        altColor 0,0,255
        parent ${samplename}_multiwig_DZNE
        
"
done


# output the individual track [container] - normalized

for i in ~/neurogen/CAGE_PDBrainMap/data_upload_scherzer_20151124/bam_files/*.plus.normalized.bw; do
    i=${i/*\//}
    samplename=${i/.plus*/}
    
echo "
    track ${samplename}_multiwignormalized_DZNE
    container multiWig
    configurable on
    shortLabel rpm_DNZE_$samplename
    longLabel DZNE normalized CAGE tags count for $samplename
    aggregate transparentOverlay
    dragAndDrop subTracks
    type bigWig 0 100
    autoScale on
    alwaysZero on
    yLineOnOff on
    yLineMark 0
    viewLimits 0:100
    visibility full
    maxHeightPixels 64:50:11
    showSubtrackColorOnUi on
    priority 1.4
    parent CAGE_BRAINCODE

        track ${samplename}_multiwignormalizedDZNE_fwd
        shortLabel DZNE normalized CAGE tags count for $samplename (fwd)
        longLabel DNZE normalized CAGE tags count for $samplename forward
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.plus.normalized.bw
        color 255,0,0
        parent ${samplename}_multiwignormalized_DZNE
        
        
        track ${samplename}_multiwignormalized_rev
        shortLabel DZNE normalized CAGE tags count for $samplename (rev)
        longLabel DZNE normalized CAGE tags count for $samplename reverse
        type bigWig 
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/cage/version3/$samplename.minus.normalized.bw
        color 0,0,255
        altColor 0,0,255
        parent ${samplename}_multiwignormalized_DZNE
"
done