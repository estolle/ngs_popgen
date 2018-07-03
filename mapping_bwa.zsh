#!/bin/zsh
############################################################################
# mapping with bwa mem and bamQC with qualimap + R viz
# Eckart Stolle, Sept 2017
############################################################################
# uses zsh shell (not bash), get it by: sudo apt-get install zsh
# usage: ./mapping_bwa.zsh.zsh INPUTFOLDER OUTPUTFOLDER CPUs REFERENCE
#fastq input: folder containing abcdefg.R1.fastq / abcdefg.R2.fastq

####uses these progz:
which bwa
which sambamba
which qualimap
which R
QUALIMAPJAR="/usr/local/src/qualimap-build-28-08-17/qualimap.jar"
ls $QUALIMAPJAR



CPUs=30
INPUTFOLDER="/scratch/scratchspace/melipona/data/merged"
OUTPUTFOLDER="/home/es/scratch/scratchspace/melipona/bam"
REF="$HOME/ref/ref_w_controls/Mqua.fa"

#cat ~/ref/Mqua/Mqua.fa ~/ref/Mqua/NC_026198_Mscut_mt.fasta ~/ref/wolbachia/wolbachia.fasta ~/ref/phiX/NC_001422.fasta > ~/ref/ref_w_controls/Mqua.fa
#bwa index ~/ref/ref_w_controls/Mqua.fa
#~/scripts/mapping_bwa.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs $REF

if [ $# -ne 4 ]; then
    echo $0: usage: ./mapping_bwa.zsh INPUTFOLDER OUTPUTFOLDER CPUs REFERENCE
	echo "\nINPUTFOLDER: Folder containing forward and reverse reads (.fastq), format aaaaa_1.fastq.gz"
	echo "\nOUTPUTFOLDER: Folder for filtered reads, will be created new if not existing"
	echo "\nCPUs: threads/cores to be used"
	echo "\nREFERENCE: reference (genome fasta) to be used"
    exit 1
fi

#conditional exit if not root
#(( $(id -u) == 0 )) || ( echo >&2 "Must be root to run script"; exit 1; )
#[[ $(id -u) == 0 ]] || ( echo >&2 "Must be root to run script"; exit 1; )
#if [[ `id -u` != 0 ]]; then
#    echo "Must be root to run script"
#exit 1
#fi

#set/get variables
INPUTFOLDER=$1
OUTPUTFOLDER=$2
CPUs=$3
REF=$4


#if some inputfiles do not exist
if [[ ! -f $(ls -1 $INPUTFOLDER/*R[1,2].fastq.gz | head -n1) ]]; then
echo "input file(s) does not exist"
    exit 1
fi
if [[ -f $(ls -1 $INPUTFOLDER/*R[1,2].fastq.gz | head -n1) ]]; then
echo "input file(s) found"
fi

if [[ ! -f $REF ]]; then
echo "reference not exist"
    exit 1
fi
if [[ -f $REF ]]; then
echo "reference found"
fi

if [[ ! -f $REF.amb ]]; then
echo "reference index not complete, creating index"
    bwa index $REF
fi


echo "mapping and bamqc"

mkdir -p $OUTPUTFOLDER
TMPFOLDER="$HOME/tmp"
mkdir -p $TMPFOLDER
mkdir -p $OUTPUTFOLDER/bamqc
#mkdir -p $OUTPUTFOLDER/bamqc_init
#mkdir -p $OUTPUTFOLDER/bamqc_dedup

#tmp dir:
#--tmpdir $HOME/overflow
mkdir -p /backup/overflow

ls -1 $INPUTFOLDER/*.R[1,2].fastq.gz | rev | cut -d"." -f 4- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/mappingsamples.lst
LISTtoPROCESS="$OUTPUTFOLDER/mappingsamples.lst"
N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
echo "processing $N samples (read-pair files)"


#dryrun
echo "settings: \n--xxxxx \nxxxxx \nxxxxx \nxxxxx \nxxxx \n-z \nxxxxxxx"
sleep 3s
cat $LISTtoPROCESS | parallel --no-notice -j 1 "echo $INPUTFOLDER/{}.R1.fastq.gz $INPUTFOLDER/{}.R2.fastq.gz"
sleep 3s

######################

N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
i=47
for (( i = 47 ; i < $N+1 ; i++))
 do
SAMPLE=$(cat $LISTtoPROCESS | sed -n $i'p')
echo $SAMPLE
##ReadgroupInfo
SAMPLEID=$SAMPLE
ID=$SAMPLE
LANE=$SAMPLE
CLADE=$SAMPLE
LIBRARY="truseq_pcrfree"



##mapping/sorting
#mkdir -p $TMPFOLDER/$SAMPLE

bwa mem -v 0 -B 6 -O 6 -E 2 -L 25 -U 50 -M -t $CPUs -R "@RG\tID:$ID\tSM:$SAMPLEID\tPL:$LANE\tLB:$LIBRARY\tPU:$CLADE" -a $REF $INPUTFOLDER/$SAMPLE.R1.fastq.gz $INPUTFOLDER/$SAMPLE.R2.fastq.gz | sambamba view --sam-input -f bam -l 0 /dev/stdin | sambamba sort -m 30G -l 9 -p -t $CPUs --tmpdir $HOME/overflow -o /dev/stdout /dev/stdin | tee $OUTPUTFOLDER/$SAMPLEID.bam | sambamba index -p -t $CPUs /dev/stdin $OUTPUTFOLDER/$SAMPLEID.bai

echo "done mapping"


#BAMQC
echo "running qualimap for bamqc"
qualimap bamqc -bam $OUTPUTFOLDER/$SAMPLEID.bam --java-mem-size=20G -nw 4000 -nr 10000 -os -hm 6 -nt $CPUs -outformat HTML -outdir $OUTPUTFOLDER/bamqc/$SAMPLEID.bamqc

echo "done QC"

#### use sambamba for dupl marking
#sambamba markdup -t $CPUs -l 8 -p --tmpdir=$TMPFOLDER --sort-buffer-size=10240 --io-buffer-size=256 --hash-table-size=1048576 --overflow-list-size=1000000 $TMPFOLDER/$SAMPLE/$SAMPLEID.bam $OUTPUTFOLDER/$SAMPLEID.bam
#sambamba index -p -t $CPUs $OUTPUTFOLDER/$SAMPLEID.bam $OUTPUTFOLDER/$SAMPLEID.bai

#qualimap bamqc -bam $OUTPUTFOLDER/$SAMPLEID.bam --java-mem-size=20G -nw 4000 -nr 10000 -os -hm 6 -nt $CPUs -outformat HTML -outdir $OUTPUTFOLDER/bamqc_dedup/$SAMPLEID.bamqc

#rm -rf $TMPFOLDER/$SAMPLE

echo $SAMPLE
done

## run bamqc-stats/viz script
$HOME/scripts/bamqc.zsh $OUTPUTFOLDER/bamqc






######### merge bams
#compress bamMerge/bamQC.tar.gz bamMerge/bamQC/

#SAMPLE="708508"
#sambamba merge -t $CPUs -l 0 -p /dev/stdout bamdedup/WTCHG_*_$SAMPLE.bam | sambamba sort -m 80G -l 5 -p -t $CPUs --tmpdir=/home/estolle/scratch/tmp3 -o /dev/stdout /dev/stdin | tee bamMerge/$SAMPLE.bam | sambamba index -p -t $CPUs /dev/stdin bamMerge/$SAMPLE.bai

#mkdir bamMerge/bamQC
#ls -1 bamMerge/*.bam | head -n30 | rev | cut -d"." -f2- | cut -d"/" -f1 | rev | parallel -j 1 "qualimap bamqc -bam bamMerge/{}.bam --java-mem-size=50G -nw 4000 -nr 10000 -os -hm 6 -nt 5 -outformat HTML -outdir bamMerge/bamQC/{}"



























