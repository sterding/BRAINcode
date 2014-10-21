# Uage:
# bash ~/neurogen/pipeline/RNAseq/modules/_make_trackDb.sh > ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.txt
# scp ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.v2.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
# to get the sample list: ls -1 ~/neurogen/rnaseq_PD/run_output/ | sed 's/_/\t/g' | sort -k3,3 -k1,1 | awk '{printf $2"="$2" "}'
#
# color selected from: http://colorbrewer2.org/

# output the general header

echo "
track RNAseq_PD_v2
shortLabel RNAseq_PD_v2
longLabel RNAseq PD (version 2)
dataVersion V2 (Oct 2014)
type bed 3
visibility hide
boxedCfg on
priority 24
compositeTrack on
dragAndDrop subTracks
subGroup1 condition Condition HC=Health_Control ILB=incidental_Lewy_body_disease PD=Parkinson's_disease AD=Alzheimer's_disease HD=Huntington's_disease
subGroup2 cellType Cell_Type SNDA=substantia_nigra_dopamine_neuron MCPY=motor_cortex_pyramidal_neuron TCPY=temporal_cortex_pyramidal_neuron
subGroup3 batch Batch batch1=1 batch2=2 batch3=3 batch4=4 batch5=5 batch6=6
subGroup4 rep Replicate rep1=1 rep2=2 rep3=3
subGroup5 view Views rawSignal=Raw_Signal rpmSignal=RPM_Signal
subGroup6 mapper Mapper uniq=Unique_Mapper multi=Multiple_Mapper
subGroup7 sample Sample NZ-H131=NZ-H131 NZ-H150=NZ-H150 NZ-H152=NZ-H152 BN00-14=BN00-14 BN03-41=BN03-41 BN04-05=BN04-05 BN04-38=BN04-38 BN05-12=BN05-12 BN06-05=BN06-05 BN07-28=BN07-28 BN07-37=BN07-37 BN08-40=BN08-40 BN08-44=BN08-44 BN08-55=BN08-55 BN08-64=BN08-64 BN08-88=BN08-88 BN09-31=BN09-31 BN09-35=BN09-35 BN10-22=BN10-22 BN10-26=BN10-26 BN10-39=BN10-39 BN10-63=BN10-63 BN10-70=BN10-70 BN10-74=BN10-74 BN11-04=BN11-04 BN11-15=BN11-15 BN11-81=BN11-81 BN11-93=BN11-93 BN11-98=BN11-98 BN12-44=BN12-44 BN97-02=BN97-02 BN97-17=BN97-17 BN99-25=BN99-25 BN99-58=BN99-58 MD4263=MD4263 MD4789=MD4789 MD5028=MD5028 MD5088=MD5088 MD5247=MD5247 MGH1000=MGH1000 MGH1000=MGH1000 MGH1026=MGH1026 NZ-H102=NZ-H102 NZ-H118=NZ-H118 NZ-H137=NZ-H137 NZ-H83=NZ-H83 UWA479=UWA479 UWA616=UWA616 BN00-34=BN00-34 BN00-53=BN00-53 BN03-50=BN03-50 BN04-64=BN04-64 BN05-33=BN05-33 BN05-61=BN05-61 BN06-25=BN06-25 BN07-26=BN07-26 BN09-05=BN09-05 BN09-34=BN09-34 BN10-10=BN10-10 BN10-90=BN10-90 BN11-05=BN11-05 BN11-28=BN11-28 BN11-41=BN11-41 BN11-60=BN11-60 BN11-86=BN11-86 BN12-25=BN12-25 BN12-28=BN12-28 BN12-54=BN12-54 BN97-53=BN97-53 BN98-32=BN98-32 BN99-50=BN99-50 BN99-54=BN99-54 BN99-66=BN99-66 MCL7798=MCL7798 MCL7842=MCL7842 MCL7878=MCL7878 BN02-04=BN02-04 BN02-24=BN02-24 BN02-33=BN02-33 BN03-15=BN03-15 BN04-52=BN04-52 BN05-10=BN05-10 BN05-10=BN05-10 BN05-16=BN05-16 BN06-57=BN06-57 BN07-11=BN07-11 BN08-85=BN08-85 BN08-90=BN08-90 BN08-90=BN08-90 BN09-20=BN09-20 BN97-10=BN97-10 BN97-37=BN97-37 BN98-19=BN98-19 BN99-44=BN99-44 BN04-42=BN04-42 BN05-17=BN05-17 BN09-12=BN09-12 BN10-28=BN10-28 BN11-54=BN11-54 BN11-90=BN11-90 BN12-33=BN12-33 BN12-42=BN12-42 BN12-55=BN12-55 BN12-56=BN12-56 BN13-05=BN13-05 BN13-17=BN13-17 BN13-18=BN13-18 MCL6444=MCL6444 MGH1288=MGH1288 MGH1288=MGH1288 MGH1488=MGH1488 MGH1488=MGH1488 NZ-PD34=NZ-PD34 NZ-PD8=NZ-PD8 UWA734=UWA734 BN03-41=BN03-41 BN05-12=BN05-12 BN06-05=BN06-05 BN08-40=BN08-40 BN10-22=BN10-22 BN10-63=BN10-63 BN10-70=BN10-70 BN11-15=BN11-15 BN12-44=BN12-44
dimensions dimX=batch dimY=sample dimA=cellType dimB=condition dimC=mapper dimD=rep
dimensionAchecked SNDA,MCPY,TCPY,
dimensionBchecked HC,PD,ILB
dimensionCchecked uniq
dimensionDchecked rep1,rep2
sortOrder condition=+ cellType=+ batch=+ rep=+ view=+ mapper=+
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
    parent RNAseq_PD_v2
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
    [ "$rep" = "" ] && rep='rep1'
    
    color="44,162,95"  # HC_SNDA
    [[ "$condition" = "PD" ]] && color="255,0,0"
    [[ "$condition" = "ILB" ]] && color="254,178,76"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "MCPY" ]] && color="153,216,201"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "TCPY" ]] && color="158,188,218"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "SNDA" ]] && color="44,162,95"
    
    
    #for j in multi uniq;
    for j in uniq;
    do
        echo "
        track ${i}_rpmSignal_${j}
        bigDataUrl http://pd:brain@panda.partners.org/~xd010/rnaseq_PD/version2/$i.$j.accepted_hits.normalized2.bw
        type bigWig
        shortLabel $i $j rpm
        longLabel RNAseq $j Normalized Signal  ( $condition $sample batch$batch $rep)
        color $color
        priority 1
        subGroups view=rpmSignal sample=$sample condition=$condition cellType=$cell batch=batch$batch rep=$rep mapper=$j
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
    parent RNAseq_PD
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