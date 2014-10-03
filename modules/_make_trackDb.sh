# Uage:
# bash ~/neurogen/pipeline/RNAseq/modules/_make_trackDb.sh > ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.txt
# scp ~/neurogen/rnaseq_PD/for_display/trackDb.RNAseq.txt xd010@panda.dipr.partners.org:~/public_html/myHub/hg19/
# to get the sample list: ls -1 ~/neurogen/rnaseq_PD/run_output/ | sed 's/_/\t/g' | sort -k3,3 -k1,1 | awk '{printf $2"="$2" "}'
#
# color selected from: http://colorbrewer2.org/

# output the general header

echo "
track RNAseq_PD
shortLabel RNAseq_PD
longLabel RNAseq PD
dataVersion Jan 2014
type bed 3
visibility hide
boxedCfg on
priority 24
compositeTrack on
dragAndDrop subTracks
subGroup1 cellType Cell_Type SNDA=substantia_nigra_dopamine_neuron MCPY=motor_cortex_pyramidal_neuron TCPY=temporal_cortex_pyramidal_neuron
subGroup2 batch Batch batch1=1 batch2=2 batch3=3 batch4=4
subGroup3 condition Condition HC=Health_Control PD=Parkinson's_disease ILB=incidental_Lewy_body_disease AD=Alzheimer's_disease HD=Huntington's_disease
subGroup4 view Views rawSignal=Raw_Signal rpmSignal=RPM_Signal
subGroup5 mapper Mapper uniq=Unique_Mapper multi=Multiple_Mapper
subGroup6 sample Sample MD5028=MD5028 MGH1000=MGH1000 MGH1026=MGH1026 MCL6444=MCL6444 MGH1288=MGH1288 MGH1488=MGH1488 BN10-22=BN10-22 BN10-39=BN10-39 MD4789=MD4789 MD5088=MD5088 MD5247=MD5247 UWA479=UWA479 UWA616=UWA616 BN09-34=BN09-34 BN97-53=BN97-53 BN98-32=BN98-32 BN99-66=BN99-66 MCL7798=MCL7798 MCL7878=MCL7878 UWA734=UWA734 BN04-05=BN04-05 BN07-28=BN07-28 BN07-37=BN07-37 BN08-44=BN08-44 BN08-55=BN08-55 BN08-88=BN08-88 BN09-31=BN09-31 BN09-35=BN09-35 BN10-26=BN10-26 BN10-63=BN10-63 BN10-70=BN10-70 BN10-74=BN10-74 BN11-04=BN11-04 BN11-15=BN11-15 BN11-81=BN11-81 BN11-93=BN11-93 BN11-98=BN11-98 BN97-17=BN97-17 MD4263=MD4263 BN03-50=BN03-50 BN04-64=BN04-64 BN05-33=BN05-33 BN05-61=BN05-61 BN09-05=BN09-05 BN10-10=BN10-10 BN11-41=BN11-41 BN11-60=BN11-60 MCL7842=MCL7842 BN00-14=BN00-14 BN03-41=BN03-41 BN04-38=BN04-38 BN05-12=BN05-12 BN06-05=BN06-05 BN08-40=BN08-40 BN08-64=BN08-64 BN12-44=BN12-44 BN97-02=BN97-02 BN99-25=BN99-25 BN99-58=BN99-58 MGH1000=MGH1000 NZ-H102=NZ-H102 NZ-H118=NZ-H118 NZ-H131=NZ-H131 NZ-H137=NZ-H137 NZ-H150=NZ-H150 NZ-H152=NZ-H152 NZ-H83=NZ-H83 BN07-26=BN07-26 BN10-90=BN10-90 BN11-86=BN11-86 BN99-50=BN99-50 BN99-54=BN99-54 MGH1288=MGH1288 MGH1488=MGH1488 NZ-PD34=NZ-PD34 NZ-PD8=NZ-PD8
dimensions dimX=batch dimY=sample dimA=cellType dimB=condition dimC=mapper
sortOrder cellType=+ condition=+ batch=+ view=+ mapper=+
dimensionAchecked SNDA
dimensionBchecked HC,PD,ILB
dimensionCchecked uniq,multi
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
    parent RNAseq_PD
"
for i in `ls -1 ~/neurogen/rnaseq_PD/run_output/`;
do
    echo "## $i";
    
    # split $i, e.g. ILB_BN03-50_SNDA_3 into 3 parts: ILB, BN03-50, SNDA, and 3
    condition=`echo $i | cut -f1 -d'_'`
    sample=`echo $i | cut -f2 -d'_'`
    cell=`echo $i | cut -f3 -d'_'`
    batch=`echo $i | cut -f4 -d'_'`
    
    color="44,162,95"  # HC_SNDA
    [[ "$condition" = "PD" ]] && color="255,0,0"
    [[ "$condition" = "ILB" ]] && color="254,178,76"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "MCPY" ]] && color="153,216,201"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "TCPY" ]] && color="229,245,249"
    [[ "$condition" = "HC" ]] && [[ "$cell" = "SNDA" ]] && color="44,162,95"
    
    
    for j in multi uniq;
    do
        echo "
        track ${i}_rpmSignal_${j}
        bigDataUrl http://panda.partners.org/~xd010/rnaseq_PD/$i.$j.accepted_hits.normalized.bw
        type bigWig
        shortLabel $i $j rpm
        longLabel RNAseq $j Normalized Signal  ( $condition $sample batch$batch )
        color $color
        priority 1
        subGroups view=rpmSignal sample=$sample condition=$condition cellType=$cell batch=batch$batch mapper=$j
        parent RNAseqRPMsignal
        "
    done
done


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