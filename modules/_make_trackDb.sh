# Uage:
# bash ~/neurogen/pipeline/RNAseq/modules/_make_trackDb.sh > ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.txt
# scp ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.v2.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
# to get the sample list: ls -1 ~/neurogen/rnaseq_PD/run_output/ | sed 's/_/\t/g' | cut -f2 | sort -u | awk '{printf $1"="$1" "}'
#
# color selected from: http://colorbrewer2.org/

# output the general header [superTrack]

echo "track RNAseq_PD
shortLabel RNA-seq
longLabel BRAINCODE RNA-seq (version 2)
dataVersion Version 2 (Oct 2014)
type bed 3
visibility full
boxedCfg on
priority 24
superTrack on show
"

# output the combined track [container multiWig]

echo "track merged_by_trimmedmean
    shortLabel combined_track
    longLabel BRAINCODE RNA-seq combined tracks by trimmed mean (5%)
    container multiWig
    aggregate none
    showSubtrackColorOnUi on
    type bigWig 0 1000
    viewLimits 0:10
    maxHeighPixels 100:50:5
    viewLimitsMax 0:30
    visibility full
    autoScale on
    alwaysZero on
    yLineOnOff on
    windowingFunction maximum
    parent RNAseq_PD
    
        track HC_MCPY_trimmedmean_uniq
        parent merged_by_trimmedmean
        shortLabel HC_MCPY trimmedmean uniq rpm
        longLabel RNAseq uniq Normalized Signal  ( HC_MCPY trimmedmean N=3)
        color 153,216,201	
        type bigWig
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/merged/HC_MCPY.trimmedmean.uniq.normalized.bw
        yLineMark 0.0238866 #7.39456e+07/3095693983
    
        track HC_TCPY_trimmedmean_uniq
        parent merged_by_trimmedmean
        shortLabel HC_TCPY trimmedmean uniq rpm
        longLabel RNAseq uniq Normalized Signal  ( HC_TCPY trimmedmean N=9)
        color 158,188,218	
        type bigWig
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/merged/HC_TCPY.trimmedmean.uniq.normalized.bw
        yLineMark 0 # 0.0238866 #7.39456e+07/3095693983  ## TO ADD
        
        track HC_SNDA_trimmedmean_uniq
        parent merged_by_trimmedmean
        shortLabel HC_SNDA trimmedmean uniq rpm
        longLabel RNAseq uniq Normalized Signal  ( HC_SNDA trimmedmean N=58)
        color 44,162,95		
        type bigWig 0 1139
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/merged/HC_SNDA.trimmedmean.uniq.normalized.bw
        yLineMark 0.01353461 #4.1899e+07/3095693983
        
        track ILB_SNDA_trimmedmean_uniq
        parent merged_by_trimmedmean
        shortLabel ILB_SNDA trimmedmean uniq rpm
        longLabel RNAseq uniq Normalized Signal  ( ILB_SNDA trimmedmean N=28)
        color 254,178,76	
        type bigWig
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/merged/ILB_SNDA.trimmedmean.uniq.normalized.bw
        yLineMark 0.0118953 #3.68243e+07/3095693983
    
        track PD_SNDA_trimmedmean_uniq
        parent merged_by_trimmedmean
        shortLabel PD_SNDA trimmedmean uniq rpm
        longLabel RNAseq uniq Normalized Signal  ( PD_SNDA trimmedmean N=21)
        color 255,0,0
        type bigWig
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/merged/PD_SNDA.trimmedmean.uniq.normalized.bw
        yLineMark 0.0179008 #5.54155e+07/3095693983
"

# output the individual track [compositeTrack]

echo "
    track RNAseq_PD_individual
    shortLabel individual_track
    longLabel BRAINCODE RNA-seq tracks for all samples
    visibility full
    boxedCfg on
    priority 24
    compositeTrack on
    dragAndDrop subTracks
    parent RNAseq_PD
    subGroup1 condition Condition HC=Health_Control ILB=incidental_Lewy_body_disease PD=Parkinson's_disease AD=Alzheimer's_disease HD=Huntington's_disease
    subGroup2 cellType Cell_Type aMCPY=motor_cortex_pyramidal_neuron aTCPY=temporal_cortex_pyramidal_neuron SNDA=substantia_nigra_dopamine_neuron
    subGroup3 batch Batch batch1=1 batch2=2 batch3=3 batch4=4 batch5=5 batch6=6
    subGroup4 rep Replicate rep1=1 rep2=2 rep3=3
    subGroup5 view Views rawSignal=Raw_Signal rpmSignal=RPM_Signal
    subGroup6 mapper Mapper uniq=Unique_Mapper multi=Multiple_Mapper
    subGroup7 sample Sample BN00-14=BN00-14 BN00-34=BN00-34 BN00-53=BN00-53 BN02-04=BN02-04 BN02-24=BN02-24 BN02-33=BN02-33 BN03-15=BN03-15 BN03-41=BN03-41 BN03-50=BN03-50 BN04-05=BN04-05 BN04-38=BN04-38 BN04-42=BN04-42 BN04-52=BN04-52 BN04-64=BN04-64 BN05-10=BN05-10 BN05-12=BN05-12 BN05-16=BN05-16 BN05-17=BN05-17 BN05-33=BN05-33 BN05-61=BN05-61 BN06-05=BN06-05 BN06-25=BN06-25 BN06-57=BN06-57 BN07-11=BN07-11 BN07-26=BN07-26 BN07-28=BN07-28 BN07-37=BN07-37 BN08-40=BN08-40 BN08-44=BN08-44 BN08-55=BN08-55 BN08-64=BN08-64 BN08-85=BN08-85 BN08-88=BN08-88 BN08-90=BN08-90 BN09-05=BN09-05 BN09-12=BN09-12 BN09-20=BN09-20 BN09-31=BN09-31 BN09-34=BN09-34 BN09-35=BN09-35 BN10-10=BN10-10 BN10-22=BN10-22 BN10-26=BN10-26 BN10-28=BN10-28 BN10-39=BN10-39 BN10-63=BN10-63 BN10-70=BN10-70 BN10-74=BN10-74 BN10-90=BN10-90 BN11-04=BN11-04 BN11-05=BN11-05 BN11-15=BN11-15 BN11-28=BN11-28 BN11-41=BN11-41 BN11-54=BN11-54 BN11-60=BN11-60 BN11-81=BN11-81 BN11-86=BN11-86 BN11-90=BN11-90 BN11-93=BN11-93 BN11-98=BN11-98 BN12-25=BN12-25 BN12-28=BN12-28 BN12-33=BN12-33 BN12-42=BN12-42 BN12-44=BN12-44 BN12-54=BN12-54 BN12-55=BN12-55 BN12-56=BN12-56 BN13-05=BN13-05 BN13-17=BN13-17 BN13-18=BN13-18 BN97-02=BN97-02 BN97-10=BN97-10 BN97-17=BN97-17 BN97-37=BN97-37 BN97-53=BN97-53 BN98-19=BN98-19 BN98-32=BN98-32 BN99-25=BN99-25 BN99-44=BN99-44 BN99-50=BN99-50 BN99-54=BN99-54 BN99-58=BN99-58 BN99-66=BN99-66 MCL6444=MCL6444 MCL7798=MCL7798 MCL7842=MCL7842 MCL7878=MCL7878 MD4263=MD4263 MD4789=MD4789 MD5028=MD5028 MD5088=MD5088 MD5247=MD5247 MGH1000=MGH1000 MGH1026=MGH1026 MGH1288=MGH1288 MGH1488=MGH1488 NZ-H102=NZ-H102 NZ-H118=NZ-H118 NZ-H131=NZ-H131 NZ-H137=NZ-H137 NZ-H150=NZ-H150 NZ-H152=NZ-H152 NZ-H83=NZ-H83 NZ-PD34=NZ-PD34 NZ-PD8=NZ-PD8 UWA479=UWA479 UWA616=UWA616 UWA734=UWA734
    dimensions dimX=batch dimY=sample dimA=cellType dimB=condition dimC=mapper dimD=rep
    dimensionAchecked SNDA,aMCPY,aTCPY,
    dimensionBchecked HC,PD,ILB
    dimensionCchecked uniq
    dimensionDchecked rep1,rep2
    sortOrder condition=+ cellType=+ batch=+ sample=+ rep=+ view=+ mapper=+
"

# output the header for normalized signal (*normalized.bw)
echo "
        track RNAseqRPMsignal
        shortLabel Normalized Signal
        view rpmSignal
        visibility full
        autoScale on
        viewLimits 0:10
        viewLimitsMax 0:30
        yLineOnOff on
        yLineMark 0
        alwaysZero on
        windowingFunction maximum
        maxHeightPixels 20:15:5
        parent RNAseq_PD_individual
"
for i in `ls -1 ~/neurogen/rnaseq_PD/run_output/`;
do
    echo "## $i";
    
    # split $i, e.g. ILB_BN03-50_SNDA_3 into 3 parts: ILB, BN03-50, SNDA, and 3 (version 1)
    # ND_BN04-52_SNDA_5c_rep1 --> ND, BN04-52, SNDA, 5c, rep1
    condition=`echo $i | cut -f1 -d'_'`
    sample=`echo $i | cut -f2 -d'_'`
    cell=`echo $i | cut -f3 -d'_'`
    batch=`echo $i | cut -f4 -d'_' | sed 's/[a-z]//'`
    rep=`echo $i | cut -f5 -d'_'`;
    
    # format
    [ "$condition" = "ND" ] && condition='HC'
    
    sampleID=$condition"_"$sample"_"$cell"_"$batch"_"$rep
    
    color="44,162,95"  # HC_SNDA
    [[ "$condition" = "PD" ]] && color="255,0,0"
    [[ "$condition" = "ILB" ]] && color="254,178,76"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "MCPY" ]] && color="153,216,201"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "TCPY" ]] && color="158,188,218"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "SNDA" ]] && color="44,162,95"
    
    cell_sorted=$cell
    [[ "$cell" = "MCPY" ]] && cell_sorted="aMCPY"
    [[ "$cell" = "TCPY" ]] && cell_sorted="aTCPY"
    
    #for j in multi uniq;
    for j in uniq;
    do
        echo "
            track ${sampleID}_rpmSignal_${j}
            bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/$i.$j.accepted_hits.normalized2.bw
            type bigWig
            shortLabel $sampleID.$j.rpm
            longLabel RNAseq $j Normalized Signal  ( $condition $sample $cell batch$batch $rep)
            color $color
            priority 1
            subGroups view=rpmSignal sample=$sample condition=$condition cellType=$cell_sorted batch=batch$batch rep=$rep mapper=$j
            parent RNAseqRPMsignal
        "
    done
done

############
exit
############

# output the header for raw signal (*bw)
echo "
        track RNAseqRawsignal
        shortLabel Raw Signal
        view rawSignal
        visibility full
        autoScale on
        viewLimits 0:10
        viewLimitsMax 0:30
        yLineOnOff on
        yLineMark 0
        alwaysZero on
        windowingFunction maximum
        maxHeightPixels 20:15:5
        parent RNAseq_PD_individual
"
for i in `ls -1 ~/neurogen/rnaseq_PD/run_output/`;
do
    echo "## $i";
    
    # split $i, e.g. ILB_BN03-50_3 into 3 parts: ILB, BN03-50 and 3
    condition=`echo $i | cut -f1 -d'_'`
    sample=`echo $i | cut -f2 -d'_'`
    cell=`echo $i | cut -f3 -d'_'`
    batch=`echo $i | cut -f4 -d'_'`
    
    color="44,162,95"  # HC_SNDA
    [[ "$condition" = "PD" ]] && color="255,0,0"
    [[ "$condition" = "ILB" ]] && color="254,178,76"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "MCPY" ]] && color="153,216,201"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "TCPY" ]] && color="153,216,201"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "SNDA" ]] && color="44,162,95"
    
    
    for j in multi uniq;
    do
        echo "
            track ${i}_RawSignal_${j}
            bigDataUrl http://panda.partners.org/~xd010/rnaseq_PD/$i.$j.accepted_hits.bw
            type bigWig
            shortLabel $i $j raw
            longLabel RNAseq $j Raw Signal  ( $condition $sample batch$batch )
            color $color
            priority 2
            subGroups view=rawSignal sample=$sample condition=$condition cellType=$cell batch=batch$batch mapper=$j
            parent RNAseqRawsignal
        "
    done
done