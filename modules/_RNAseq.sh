###########################################
# bash script for running paired-end RNAseq
# author: Xianjun Dong
# email: xdong@rics.bwh.harvard.edu
# date: 9/16/2013
# version: 2.0
# Note: call this script in the folder of fastq file
## TODO:
#1. use PINSEQ for QC
#2. use STAR for mapping
###########################################
#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: $0 HC_BN10-39_2.R1.fastq.gz HC_BN10-39_2.R2.fastq.gz"
  exit
fi

###########################################
echo "["`date`"] STEP 1. Configuring"
###########################################

modulename=`basename $0`
set +o posix  #  enables the execution of process substitution e.g. http://www.linuxjournal.com/content/shell-process-redirection

R1=$1  # filename of R1 
R2=$2  # filename of R2 (for paired-end reads)
samplename=${R1/[.|_]R1*/}

pipeline_path=$HOME/neurogen/pipeline/RNAseq
source $pipeline_path/config.txt

strandoption="--library-type fr-unstranded"; split="-nosplit"; 
[[ $samplename == *stranded* ]] && strandoption="--library-type fr-secondstrand"  # by default we use Illumina SMARTer stranded RNA-Seq kit
[[ $samplename == *stranded* ]] && split="-split"

inputdir=$PWD
outputdir=$inputdir/../run_output
[ -d $outputdiri/$samplename ] || mkdir -p $outputdir/$samplename

###########################################
echo "["`date`"] STEP 2. quality filter: adaptor removal/clip"
###########################################
## TODO: Reads with PolyA/T tail were not trimmed. Using PRINSEQ instead for removal

##### adaptor removal
[ -d $inputdir/../filtered ] || mkdir $inputdir/../filtered

[ -e $inputdir/../filtered/adaptor.fa ] || echo -e ">Truseq_complementary_part_3p\nAGATCGGAAGAGC" > $inputdir/../filtered/adaptor.fa

# remove adaptor with fastq-mcf (https://code.google.com/archive/p/ea-utils/wikis/FastqMcf.wiki)
[ ! -f $outputdir/$samplename/.status.$modulename.adaptorremoval ] && \
fastq-mcf -o $inputdir/../filtered/$R1 -o $inputdir/../filtered/$R2 -t 0 -x 10 -l 15 -w 4 -q 10 -u $inputdir/../filtered/adaptor.fa <(zcat $R1) <(zcat $R2) && \
touch $outputdir/$samplename/.status.$modulename.adaptorremoval 

cd $inputdir/../filtered

#############################################
echo "["`date`"] STEP 2.1. extract short reads (e.g. <=37nt, R1/R2 same length and reverse complementary)"
############################################
[ -d $inputdir/../filtered/shortReads ] || mkdir $inputdir/../filtered/shortReads

[ ! -f $outputdir/$samplename/.status.$modulename.shortReadsExtract ] && \
#zcat $R1 | perl -ne '$h=$_; $s=<>; $a=<>; $ss=<>; print "$h$s$a$ss" if(length($s)<=31);' | gzip > shortReads/$R1 
#zcat $R2 | perl -ne '$h=$_; $s=<>; $a=<>; $ss=<>; print "$h$s$a$ss" if(length($s)<=31);' | gzip > shortReads/$R2 
join -j 1 <(zcat $R1 | perl -ne 'chomp;$h=$_;chomp($s=<>);chomp($a=<>);$ss=<>; print "$h\t$s\t$a\t$ss" if(length($s)<=37);' | sort -k1,1) <(zcat $R2 | perl -ne 'chomp;$h=$_;chomp($s=<>);chomp($a=<>);$ss=<>; print "$h\t$s\t$a\t$ss" if(length($s)<=37);' | sort -k1,1) | perl -ne '@a=split;$rcs=reverse($a[2]);$rcs=~tr/ACGTacgt/TGCAtgca/; print ">$a[0]\n$a[2]\n" if($rcs eq $a[6]);' | gzip > shortReads/${R1/R1.fastq/fasta} && \
touch $outputdir/$samplename/.status.$modulename.shortReadsExtract 

[ ! -f $outputdir/$samplename/.status.$modulename.shortReadsExtract_stemloop ] && \
zcat shortReads/${R1/R1.fastq/fasta} | perl -ne '$h=$_; $s=<>; chomp($s); $sub=reverse(substr($s,-6)); $sub=~tr/ACGTacgt/TGCAtgca/; print "$h$s\n" if($s!~/GCGCGCGCGC/ && (substr($s,0,6) eq $sub));' | gzip > shortReads/${R1/R1.fastq/fasta.stemloop} && \
touch $outputdir/$samplename/.status.$modulename.shortReadsExtract_stemloop

#############################################
echo "["`date`"] STEP 3. QC"
############################################

[ ! -f $outputdir/$samplename/.status.$modulename.fastqc ] && \
fastqc --outdir $outputdir/$samplename --extract -t 2 $R1 $R2 && \
rm $outputdir/$samplename/*fastqc.zip && \
touch $outputdir/$samplename/.status.$modulename.fastqc

#############################################
echo "["`date`"] STEP 3.1 k-mer for outlier detection"
############################################
# require to install kpal: http://kpal.readthedocs.org/en/latest/install.html
[ ! -f $outputdir/$samplename/.status.$modulename.kmers ] && \
fqHEAD=`zcat $R1 | head -n1 | sed 's/@//g'| cut -f1 -d":"` && \
rm $outputdir/$samplename/*k9 && \
fastqToFa -nameVerify=$fqHEAD $R1 stdout | kpal count -p ${samplename//./_} -k 9 - $outputdir/$samplename/R1.k9 && \
fastqToFa -nameVerify=$fqHEAD $R2 stdout | kpal count -p ${samplename//./_} -k 9 - $outputdir/$samplename/R2.k9 && \
kpal merge $outputdir/$samplename/R1.k9 $outputdir/$samplename/R2.k9 $outputdir/$samplename/k9 && \
rm $outputdir/$samplename/R1.k9 $outputdir/$samplename/R2.k9 && \
touch $outputdir/$samplename/.status.$modulename.kmers

############################################
echo "["`date`"] STEP 4. mapping"
############################################
cd $inputdir/../filtered
## tophat (TODO: switch to STAR)
## Note: 
## 1) output sorted accepted_hits.bam, allow up to 100 multiple hits
## 2) GATK requires @RG group has fields of ID, SM, PL, LB, PU

## NOTE: this was done for all batch 1-5 samples with tophat version 2.0.8
[ ! -f $outputdir/$samplename/.status.$modulename.mapping ] && \
tophat -o $outputdir/$samplename --rg-id $samplename --rg-sample $samplename --rg-platform ILLUMINA --rg-library $samplename --rg-platform-unit $samplename --keep-fasta-order -p $CPU --read-mismatches $mm $tophat $PE_option $strandoption --max-multihits $MAX_HIT --no-coverage-search $GENOME/Sequence/Bowtie2Index/genome $R1 $R2 && \
touch $outputdir/$samplename/.status.$modulename.mapping

## Using BWA-MEM to prepare alignment for calling cirRNA with CIRI
#bwa mem $R1 $R2

############################################
echo "["`date`"] STEP 4.1. mapping & calling circRNA"
############################################
cd $outputdir/$samplename

mindist=200

# circular RNA calling and quantification (see: https://github.com/YangLab/CIRCexplorer/blob/master/README.md)
## using tophat-2.0.10 as suggested by CIRCexplorer
## the following option is added by Rebeca:
## --fusion-min-dist default is 10,000,000 changed to --fusion-min-dist 200
## --fusion-ignore-chromosomes added mitochondria
[ ! -f .status.$modulename.circRNA ] && \
bamToFastq -i unmapped.bam -fq unmapped.fastq &&
tophat -o tophat_fusion -p $CPU --fusion-search --fusion-min-dist $mindist --fusion-ignore-chromosomes chrM --keep-fasta-order --bowtie1 --no-coverage-search $BOWTIE_INDEXES/genome unmapped.fastq &&
CIRCexplorer.py -f tophat_fusion/accepted_hits.bam -g $GENOME/Sequence/WholeGenomeFasta/genome.fa -r $GENOME/Annotation/Genes/refFlat.txt -o circ.txt && \
touch .status.$modulename.circRNA


## update to CIRCExplore2 to include CDR1as (see email with Xiao-Ou Zhou) - 10/23/2017
## extract CDR1as gtf from UCSC and converted to genePred then manually added to refFlat.txt
[ ! -f .status.$modulename.circRNA.circexplore2 ] && \
CIRCexplorer2 parse -t TopHat-Fusion tophat_fusion/accepted_hits.bam -b back_spliced_junction.bed > .CIRCexplorer2_parse.log && \
CIRCexplorer2 annotate -r $GENOME/Annotation/Genes/refFlat_plus_CDR1as.txt -g $GENOME/Sequence/WholeGenomeFasta/genome.fa -b back_spliced_junction.bed -o circularRNA_known.txt --low-confidence > .CIRCexplorer2_annotate.log && \
touch .status.$modulename.circRNA.circexplore2

# ## calling cirRNA using Mapsplice
# [ ! -f $outputdir/$samplename/.status.$modulename.circRNA_Mapsplice ] && \
# mkdir $outputdir/$samplename/mapsplice_out && \
# python /PHShome/xd010/MapSplice-v2.1.8/mapsplice.py -p $CPU -o $outputdir/$samplename/mapsplice_out --bam --non-canonical --fusion-non-canonical --min-fusion-distance $mindist --gene-gtf $ANNOTATION_GTF -c $GENOME/Sequence/Chromosomes -x $BOWTIE_INDEXES/genome -1 $outputdir/$samplename/unmapped.fastq  && \
# rm $outputdir/$samplename/unmapped.fastq && \
# touch $outputdir/$samplename/.status.$modulename.circRNA_Mapsplice
# 
# #-p $CPU is threads
# #-c <Reference_Sequence folder>, each file in 1 chr in fasta
# #-x <Bowtie_Index folder> each file in 1 chr in fasta
# #--min-fusion-distance Minimim distance between two segments to be considered as fusion candidate. Manual recommends 200 for circRNAs
# #--min-map-len <int> leaving at default 50
# # --gene-gtf <string> Gene annotation file in GTF format, used to annotate fusion junctions. Required for the detection of Circular RNA
# # --fusion | --fusion-non-canonical for circRNAs
# 


###########################################
echo "["`date`"] STEP 5. post-processing, format converting"
###########################################

cd $outputdir/$samplename

## NOTE: after changing to output sorted bam from Tophat directly, the following sam-->bam and sorting is not needed any more, only index is retained.
[ ! -f .status.$modulename.sam2bam ] && \
#samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted && \
#mv accepted_hits.sorted.bam accepted_hits.bam && \
samtools index accepted_hits.bam && \
touch .status.$modulename.sam2bam

## ONLY NEEDED FOR GTEx RNAseq data (which has different chromosome names)
#[ ! -f .status.$modulename.rehead ] && \
#_header_change.sh && \
#touch .status.$modulename.rehead

[ ! -f .status.$modulename.bam2stat ] && \
echo `samtools view -cF 0x100 accepted_hits.bam` "primary alignments (from samtools view -cF 0x100)" > accepted_hits.bam.stat && \
samtools flagstat accepted_hits.bam >> accepted_hits.bam.stat && \
touch .status.$modulename.bam2stat

[ ! -f $outputdir/$samplename/.status.$modulename.bam2annotation ] && \
_bam2annotation.sh accepted_hits.bam > accepted_hits.bam.bam2annotation && \
Rscript $pipeline_path/modules/_bam2annotation.r accepted_hits.bam.bam2annotation accepted_hits.bam.bam2annotation.pdf && \
touch $outputdir/$samplename/.status.$modulename.bam2annotation

### bigwig for UCSC
echo "## generating bigwig files for UCSC display"

## NEW VERSION (rerun for BATCH 1-7)
[ ! -f $outputdir/$samplename/.status.$modulename.sam2bw ] && \
echo "## normalization using total reads mapped to non_rRNA_mt" && \
total_mapped_reads=`grep -w total_non_rRNA_mt accepted_hits.bam.bam2annotation | cut -f2 -d' '` && \
bam2bigwig.sh accepted_hits.bam $split $total_mapped_reads && \
#rm accepted_hits.bed? !(*normalized*).bedGraph && \
touch $outputdir/$samplename/.status.$modulename.sam2bw

# % of genome coverage vs. RNAseq density (e.g. >1rpm)
[ ! -f .status.$modulename.rpm_vs_coverage ] && \
awk 'BEGIN{max=100; UNIT=0.01; OFS="\t";}{if($0~/^#/) {print; next;} i=int($4/UNIT);if(i>max) i=max; rpm[i]+=($3-$2);}END{for(x=max;x>=0;x--) print x*UNIT, rpm[x]?rpm[x]:0;}' accepted_hits.normalized.bedGraph > accepted_hits.normalized.rpm_vs_coverage.txt && \
touch .status.$modulename.rpm_vs_coverage

#rm accepted_hits.bed accepted_hits.*bedGraph

###########################################
echo "["`date`"] STEP 6. call variation [OFF]"
###########################################

cd $outputdir/$samplename
#
#[ ! -f .status.$modulename.callSNP ] && \
#_callSNP.sh accepted_hits.sam && \
#touch .status.$modulename.callSNP

### shuilin's GATK
##[ ! -f .status.$modulename.callSNP_GATK ] && \
##_bam2vcf.sh accepted_hits.sam && \
##touch .status.$modulename.callSNP_GATK

###########################################
echo "["`date`"] STEP 7. assembly and quantification"
###########################################

cd $outputdir/$samplename

##echo "## run trinity to do de-novo discovery"
#Trinity.pl --output denovo --seqType fq --JM 100G --left $R1 --right $R2 --CPU $CPU

echo "## run cufflinks for do de-novo discovery" 
[ ! -f .status.$modulename.cufflinks.multi.denovo ] && \
cufflinks --no-update-check --no-faux-reads $strandoption -o ./denovo -p $CPU -g $ANNOTATION_GTF -M $MASK_GTF accepted_hits.bam 2> cufflinks.denovo.log && \
touch .status.$modulename.cufflinks.multi.denovo

#echo "## run Scripture to do de-novo discovery" [TURN OFF due to the possible memory leak]
## see demo at: http://garberlab.umassmed.edu/data/RNASeqCourse/analysis.workshop.pdf
## "Scripture initially transforms the genome into a graph topology, which represents all possible connections of bases in the transcriptome either when they occur consecutively or when they are connected by a spliced read. Scripture uses this graph topology to reduce the transcript reconstruction problem to a statistical segmentation problem of identifying significant transcript paths across the graph" (http://www.nature.com.ezp-prod1.hul.harvard.edu/nmeth/journal/v8/n6/full/nmeth.1613.html)
#[ ! -f .status.$modulename.scripture ] && \
#> .scripture.paraFile && \
#for i in `seq 1 22` X Y M; do echo "java -Xmx40000m -jar $pipeline_path/bin/scripture-beta2.jar -alignment accepted_hits.bam -out scripture.segments.$i -sizeFile $GENOME/Sequence/WholeGenomeFasta/hg19.chrom.size -chr chr$i -chrSequence $GENOME/Sequence/Chromosomes/chr$i.fa &>>scripture.log" >> .scripture.paraFile; done && \
#rm -f .scripture.paraFile.completed && \
#ParaFly -c .scripture.paraFile -CPU 8 && \
#touch .status.$modulename.scripture

echo "## run cufflinks to get FPKM"
# Using gtf from deno assembly
# Note: "-b" option (for bias correction) can lead to segementation fault error.
[ ! -f .status.$modulename.cufflinks ] && \
cufflinks --no-update-check $strandoption -o ./ -p $CPU -G $ANNOTATION_GTF -M $MASK_GTF --compatible-hits-norm --multi-read-correct accepted_hits.bam && \
touch .status.$modulename.cufflinks 

#echo "## run cufflinks without -M option"
#cufflinks -q --no-update-check $strandoption -o ./cufflink_wo_M -p $CPU -G $ANNOTATION_GTF -b $BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam

echo "## run htseq for reads count"
[ ! -f .status.$modulename.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q -f bam accepted_hits.bam $ANNOTATION_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch .status.$modulename.htseqcount

echo "## quantification for meta exons"
# script to generate meta exons (see README in $ANNOTATION for how to generate meta exons)
# awk '{OFS="\t"; split($4,a,"___"); split(a[4],b,"."); $4=a[1]"___"a[2]"___"b[1]; print}' exons.bed | sort -k4,4 -k1,1 -k2,2n -k3,3n | awk '{OFS="\t"; if(id!=$4 || $2>e) {if(id!="") print chr,s,e,id,1,str; chr=$1;s=$2;e=$3;id=$4;str=$6;} else if($3>e) e=$3;}END{print chr,s,e,id,1,str;}' > exons.meta.bed
[ ! -f .status.$modulename.metaexon ] && \
awk '{OFS="\t"; $4=$4"__"$1"_"$2"_"$3; print}' $ANNOTATION/exons.meta.bed | /source/bedtools2/2.18.2/coverageBed -abam accepted_hits.bam -b - -s -counts > readscount.by.metaexon.tab 2> readscount.by.metaexon.tab.stderr && \
touch .status.$modulename.metaexon

############################################
echo "["`date`"] STEP 8. for uniq"
############################################

## extract the unique mapper only
[ -d $outputdir/$samplename/uniq ] || mkdir $outputdir/$samplename/uniq
cd $outputdir/$samplename/uniq


[ ! -f $outputdir/$samplename/.status.$modulename.uniq ] && \
samtools view -H $outputdir/$samplename/accepted_hits.bam > accepted_hits.sam && \
samtools view $outputdir/$samplename/accepted_hits.bam | fgrep -w NH:i:1 >> accepted_hits.sam && \
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted && \
rm accepted_hits.sam && \
mv accepted_hits.sorted.bam accepted_hits.bam && \
samtools index accepted_hits.bam && \
cufflinks --no-update-check $strandoption -o ./ -p $CPU -G $ANNOTATION_GTF -M $MASK_GTF --compatible-hits-norm accepted_hits.bam && \
touch $outputdir/$samplename/.status.$modulename.uniq

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.bam2stat ] && \
echo `samtools view -cF 0x100 accepted_hits.bam` "primary alignments (from samtools view -cF 0x100)" > accepted_hits.bam.stat && \
samtools flagstat accepted_hits.bam >> accepted_hits.bam.stat && \
touch $outputdir/$samplename/.status.$modulename.uniq.bam2stat

echo "## run cufflinks for do de-novo discovery using uniq mapper only" 
[ ! -f $outputdir/$samplename/.status.$modulename.cufflinks.denovo ] && \
cufflinks --no-update-check --no-faux-reads $strandoption -o ./denovo -p $CPU -g $ANNOTATION_GTF -M $MASK_GTF accepted_hits.bam 2> cufflinks.denovo.log && \
touch $outputdir/$samplename/.status.$modulename.cufflinks.denovo

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.bam2annotation ] && \
_bam2annotation.sh accepted_hits.bam > accepted_hits.bam.bam2annotation && \
Rscript $pipeline_path/modules/_bam2annotation.r accepted_hits.bam.bam2annotation accepted_hits.bam.bam2annotation.pdf && \
touch $outputdir/$samplename/.status.$modulename.uniq.bam2annotation

[ -f accepted_hits.bam.bam2annotation ] || _bam2annotation.sh accepted_hits.bam > accepted_hits.bam.bam2annotation
[ ! -f $outputdir/$samplename/.status.$modulename.uniq.sam2bw ] && \
echo "## normalizing: instead of using total reads, use reads only mapped to non-rRNA-mtRNA for normalization" && \
total_mapped_reads2=`grep -w total_non_rRNA_mt accepted_hits.bam.bam2annotation | cut -f2 -d' '` && \
bam2bigwig.sh accepted_hits.bam $split $total_mapped_reads2 && \
rm accepted_hits.bed? accepted_hits*s.bedGraph && \
touch $outputdir/$samplename/.status.$modulename.uniq.sam2bw

## DEPRECATED: Now all normalized.bw are changed to use total reads mapped to non-rRNA_mt 
# [ ! -f $outputdir/$samplename/.status.$modulename.uniq.normalize ] && \
# total_mapped_reads2=`grep total_non_rRNA_mt accepted_hits.bam.bam2annotation | cut -f2 -d' '` && \
# awk -v tmr=$total_mapped_reads2 'BEGIN{OFS="\t"; print "# total-rRNA-chrM="tmr;}{$4=$4*1e6/tmr; print}' accepted_hits.bedGraph > accepted_hits.normalized2.bedGraph && \
# bedGraphToBigWig accepted_hits.normalized2.bedGraph $ANNOTATION/ChromInfo.txt accepted_hits.normalized2.bw && \
# touch $outputdir/$samplename/.status.$modulename.uniq.normalize

echo "## calcualte RPKM (based on TOTAL reads mapped to nuclear genome)"
[ ! -f $outputdir/$samplename/.status.$modulename.uniq.cufflinks.rpkm ] && \
cufflinks --no-update-check $strandoption -o ./rpkm -p $CPU -G $ANNOTATION_GTF accepted_hits.bam.non-rRNA-mt.bam 2> cufflinks.rpkm.log && \
touch $outputdir/$samplename/.status.$modulename.uniq.cufflinks.rpkm

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.cuffquant.rpkm ] && \
cuffquant --no-update-check $strandoption -o ./rpkm -p $CPU $ANNOTATION_GTF accepted_hits.bam.non-rRNA-mt.bam 2> cuffquant.rpkm.log && \
touch $outputdir/$samplename/.status.$modulename.uniq.cuffquant.rpkm

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.htseqcount ] && \
htseq-count -m intersection-strict -t exon -i gene_id -s no -q -f bam accepted_hits.bam $ANNOTATION_GTF > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
touch $outputdir/$samplename/.status.$modulename.uniq.htseqcount

[ ! -f $outputdir/$samplename/.status.$modulename.uniq.metaexon ] && \
coverageBed -abam accepted_hits.bam.non-rRNA-mt.bam -b $ANNOTATION/exons.meta.bed -s -counts > readscount.by.metaexon.tab 2> readscount.by.metaexon.tab.stderr && \
touch $outputdir/$samplename/.status.$modulename.uniq.metaexon

echo "## quantification for meta intron"
# see README in $ANNOTATION for how to generate meta introns
[ ! -f $outputdir/$samplename/.status.$modulename.uniq.metaintron ] && \
awk '{OFS="\t"; $4=$4"="$1"_"$2"_"$3; print}' $ANNOTATION/introns.meta.bed | bigWigAverageOverBed accepted_hits.normalized.bw stdin stdout 2> /dev/null | sed 's/=[^\t]*//g' | groupBy -g 1 -c 2,4 -o sum,sum | sed 's/___/\t/g' | cut -f2,4- | awk 'BEGIN{OFS="\t";print "id","intronLength","totalRPM","meanRPM"}{print $0,$3/$2}' > meanRPM.of.metaintron.by.gene.tab && \
#awk '{OFS="\t"; $4=$4"="$1"_"$2"_"$3; print}' $ANNOTATION/introns.meta.bed | coverageBed -a - -b accepted_hits.bam.non-rRNA-mt.bam -s -counts > readscount.by.metaintron.tab 2> readscount.by.metaintron.tab.stderr && \
touch $outputdir/$samplename/.status.$modulename.uniq.metaintron

#[ ! -f $outputdir/$samplename/.status.$modulename.uniq.callSNP ] && \
#_callSNP.sh accepted_hits.sam && \
#touch $outputdir/$samplename/.status.$modulename.uniq.callSNP

## shuilin's GATK
#[ ! -f $outputdir/$samplename/.status.$modulename.uniq.callSNP_GATK ] && \
#samtools view -SH $outputdir/$samplename/accepted_hits.sam > accepted_hits.sam && \
#fgrep -w NH:i:1 $outputdir/$samplename/accepted_hits.sam >> accepted_hits.sam && \
#_bam2vcf.sh accepted_hits.sam && \
#touch $outputdir/$samplename/.status.$modulename.uniq.callSNP_GATK

#[ ! -f $outputdir/$samplename/.status.$modulename.uniq.rpm_vs_coverage ] && \
#awk 'BEGIN{max=100; UNIT=0.01; OFS="\t";}{if($0~/^#/) {print; next;} i=int($4/UNIT);if(i>max) i=max; rpm[i]+=($3-$2);}END{for(x=max;x>=0;x--) print x*UNIT, rpm[x]?rpm[x]:0;}' accepted_hits.normalized.bedGraph > accepted_hits.normalized.rpm_vs_coverage.txt && \
#touch $outputdir/$samplename/.status.$modulename.uniq.rpm_vs_coverage

#rm accepted_hits.bed accepted_hits.*bedGraph


############################################
echo "["`date`"] STEP 9. prepare for tracks files to display on UCSC / IGV"
############################################
#
[ -d $inputdir/../for_display ] || mkdir $inputdir/../for_display
cd $inputdir/../for_display

### others
[ ! -f $outputdir/$samplename/.status.$modulename.makelinks ] && \
# make soft link
# multi
for i in $outputdir/$samplename/{accepted_hits.bam,accepted_hits.bam.bai,*tracking,hgseqcount.by.gene.tab,transcripts.gtf,*.bw}; do ii=${i/.*\//}; ln -fs $i $samplename.multi.$ii; done && \

# uniq
for i in $outputdir/$samplename/uniq/{accepted_hits.bam,accepted_hits.bam.bai,*tracking,hgseqcount.by.gene.tab,transcripts.gtf,*.bw}; do ii=${i/.*\//}; ln -fs $i $samplename.uniq.$ii; done && \

## QC
ln -fs $outputdir/$samplename/$samplename.R1_fastqc $samplename.R1_fastqc && \
ln -fs $outputdir/$samplename/$samplename.R2_fastqc $samplename.R2_fastqc && \

touch $outputdir/$samplename/.status.$modulename.makelinks

echo "["`date`"] DONE: $modulename job for sample $samplename is done !!"