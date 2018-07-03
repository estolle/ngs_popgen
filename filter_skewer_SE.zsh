#!/bin/zsh
############################################################################
# readfiltering with skewer
# Eckart Stolle, Sept 2017
############################################################################
# uses zsh shell (not bash), get it by: sudo apt-get install zsh
# usage: ./skewer.zsh INPUTFOLDER OUTPUTFOLDER CPUs

# !!! comment/uncomment parts if sickle is used for further/stricter trimming in addition  (trimming/stats/plot section)
# !!! see below and update location of fwd/rev adapters to filter
# !!! see below to adapt filter_QCplot.R file location
# !!! see below to adapt the format of the fwd/rev reads file formats (.R1 or _1)
# !!! L=80 length threshold
# !!! optical distance threshold depedning on platform

FWD_ADAPTERS="/scratch/genomes/illumina.fwd.adapters.fa"
REV_ADAPTERS="/scratch/genomes/illumina.rev.adapters.fa"
#AVG_QUAL_THRESHOLD=30
#END_QUAL_THRESHOLD=20
#L=80

AVG_QUAL_THRESHOLD=30
END_QUAL_THRESHOLD=20
L=80
OPTICALDIST=40
#PLOTSCRIPT="$HOME/scripts/filter_QCplot.R"
PLOTSCRIPT="$HOME/scripts/filter_QCplot_no-sickle_no-clumpify.R"



# pip install multiqc
# brew install fastqc
# install skewer, sickle, bbmap(clumpify.sh)
# for QC plots: R, packages: tidyverse (ggplot2,readr,dplyr), cowplot,scales

#./skewer.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs

#INPUTFOLDER="/scratch/ek/reads/apis/Cridland2017"
#OUTPUTFOLDER="/scratch/ek/reads/apis/Cridland2017/filtered_fastq"
#CPUs=30


#run:
#~/scripts/skewer.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs


if [ $# -ne 3 ]; then
    echo $0: usage: ./skewer.zsh INPUTFOLDER OUTPUTFOLDER CPUs 
	echo "\nINPUTFOLDER: Folder containing forward and reverse reads (.fastq), format aaaaa_1.fastq.gz"
	echo "\nOUTPUTFOLDER: Folder for filtered reads, will be created new if not existing"
	echo "\nCPUs: threads/cores to be used"
    exit 1
fi

#set/get variables
INPUTFOLDER=$1
OUTPUTFOLDER=$2
CPUs=$3

#if adapterfiles do not exist
if [[ ! -f $FWD_ADAPTERS ]]; then
echo "FWD adapter file does not exist"
fi
if [[ ! -f $REV_ADAPTERS ]]; then
echo "FWD adapter file does not exist"
fi

echo "filtering reads (skewer) and fastqc"
which skewer


mkdir -p $OUTPUTFOLDER
ls -1 $INPUTFOLDER/*_1.fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/samples.lst
LISTtoPROCESS="$OUTPUTFOLDER/samples.lst"
N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
echo "processing $N samples (read-pair files)"

#filtering
#dryrun
if [[ -f $FWD_ADAPTERS ]]; then
echo "FWD adapter file found"
else echo "FWD adapter file not found"
fi
if [[ -f $REV_ADAPTERS ]]; then
echo "REV adapter file found"
else echo "REV adapter file not found"
fi
sleep 3s
cat $LISTtoPROCESS | parallel --no-notice -j 1 "echo $INPUTFOLDER/{}'_1'.fastq.gz"


#######
TMPFOLDER=$HOME/tmp
mkdir -p $TMPFOLDER

mkdir -p $OUTPUTFOLDER/QC
mkdir -p $OUTPUTFOLDER/fastqcs
mkdir -p $OUTPUTFOLDER/multiqc

#mkdir -p $OUTPUTFOLDER/strict
#mkdir -p $OUTPUTFOLDER/strict/fastqcs
#mkdir -p $OUTPUTFOLDER/strict/multiqc

mkdir -p $INPUTFOLDER/QC
mkdir -p $INPUTFOLDER/fastqcs
mkdir -p $INPUTFOLDER/multiqc


echo "___________________"


echo "skewer/sickle parameters: \n--mean-quality $AVG_QUAL_THRESHOLD \n--end-quality $END_QUAL_THRESHOLD \n-l $L bp"
sleep 3s
echo "starting filtering"

#dupedist instead dist HiSeq2000:40 HiSeq4000:2500
#main filtering (clumpify, skewer)

#clumpify.sh in1=$INPUTFOLDER/{}'_1'.fastq.gz out1=$TMPFOLDER/{}.temp.R1.fastq.gz dedupe optical dupedist=$OPTICALDIST -Xmx100g passes=1 subs=1 k=31 spantiles=f 2> >(tee $OUTPUTFOLDER/{}.optical_duplicates.log >&1); echo clumpify_done;\

cat $LISTtoPROCESS | parallel --no-notice -j 3 "echo {}; skewer -m any $INPUTFOLDER/{}'_1'.fastq.gz -x $FWD_ADAPTERS -y $REV_ADAPTERS --mean-quality $AVG_QUAL_THRESHOLD --end-quality $END_QUAL_THRESHOLD -l $L -n -r 0.1 -z --format auto -t $CPUs -o $TMPFOLDER/{} ; echo skewer_done;\
mv -f $TMPFOLDER/{}'-trimmed'.fastq.gz $OUTPUTFOLDER/{}.R1.fastq.gz; mv -f $TMPFOLDER/{}'-trimmed'.log $OUTPUTFOLDER/{}.skewer.log ; echo {}"

#strict filtering (sickle) and moving of files #adapt/test if used for SE
#cat $LISTtoPROCESS | parallel --no-notice -j 15 "echo {}; sickle pe -t sanger --gzip-output -q $AVG_QUAL_THRESHOLD -l $L -f $TMPFOLDER/{}.R1.fastq.gz -r $TMPFOLDER/{}.R2.fastq.gz -o $OUTPUTFOLDER/strict/{}.R1.fastq.gz -p $OUTPUTFOLDER/strict/{}.R2.fastq.gz -s $OUTPUTFOLDER/strict/{}.singletons.fastq.gz 1> >(tee $OUTPUTFOLDER/strict/{}.sickle.log >&1) ; echo sickle_done ; \
#mv -f $TMPFOLDER/{}.R1.fastq.gz $OUTPUTFOLDER/; mv -f $TMPFOLDER/{}.R2.fastq.gz $OUTPUTFOLDER/; mv -f $TMPFOLDER/{}.skewer.log $OUTPUTFOLDER/; \

#rm -f $TMPFOLDER/{}.temp.R1.fastq.gz; echo {}"

### run fastqc
#fastqc --noextract -f fastq --threads $CPUs --nogroup -o $OUTPUTFOLDER/fastqcs/ $TMPFOLDER/{}.R1.fastq.gz ; fastqc --noextract -f fastq --threads $CPUs --nogroup -o $OUTPUTFOLDER/fastqcs/ $TMPFOLDER/{}.R2.fastq.gz

ls -1 $INPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel --no-notice -j 10 "echo {}; fastqc --noextract -f fastq --threads 15 --nogroup -o $INPUTFOLDER/fastqcs/ $INPUTFOLDER/{}.fastq.gz; echo {}"

ls -1 $OUTPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel --no-notice -j 10 "echo {}; fastqc --noextract -f fastq --threads 15 --nogroup -o $OUTPUTFOLDER/fastqcs/ $OUTPUTFOLDER/{}.fastq.gz; echo {}"

#ls -1 $OUTPUTFOLDER/strict/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel --no-notice -j 15 "echo {}; fastqc --noextract -f fastq --threads 5 --nogroup -o $OUTPUTFOLDER/strict/fastqcs/ $OUTPUTFOLDER/strict/{}.fastq.gz; echo {}"

### summarize fastqc's with multiqc
multiqc -o $INPUTFOLDER/multiqc/ $INPUTFOLDER/fastqcs/
multiqc -o $OUTPUTFOLDER/multiqc/ $OUTPUTFOLDER/fastqcs/
#multiqc -o $OUTPUTFOLDER/strict/multiqc/ $OUTPUTFOLDER/strict/fastqcs/

#compress
#tar -czf $INPUTFOLDER/fastqcs.tar.gz $INPUTFOLDER/fastqcs/
#tar -czf $OUTPUTFOLDER/fastqcs.tar.gz $OUTPUTFOLDER/fastqcs/
#tar -czf $OUTPUTFOLDER/strict/fastqcs.tar.gz $OUTPUTFOLDER/strict/fastqcs/

tar -cv $INPUTFOLDER/fastqcs/ | pigz -p $CPUs > $INPUTFOLDER/fastqcs.tar.gz
tar -cv $OUTPUTFOLDER/fastqcs/ | pigz -p $CPUs > $OUTPUTFOLDER/fastqcs.tar.gz
#tar -cv $OUTPUTFOLDER/strict/fastqcs/ | pigz -p $CPUs > $OUTPUTFOLDER/strict/fastqcs.tar.gz

mv -f $INPUTFOLDER/fastqcs.tar.gz $INPUTFOLDER/QC
mv -f $INPUTFOLDER/multiqc/ $INPUTFOLDER/QC

mv -f $OUTPUTFOLDER/fastqcs.tar.gz $OUTPUTFOLDER/QC
mv -f $OUTPUTFOLDER/multiqc/ $OUTPUTFOLDER/QC


rm -rf $OUTPUTFOLDER/fastqcs/ $OUTPUTFOLDER/strict/fastqcs/ $INPUTFOLDER/fastqcs/


###################### STATS
##summarize stats
N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
###



echo "SAMPLE\tREADSIN\tREADPAIRSbefore\tREADPAIRSafter\tREADPAIRSafterPC\tREADPAIRSafterTRIMMED\tREADPAIRSafterUNTRIMMED\tSICKLEOUTPAIRED\tSICKLEPC\tSICKLEOUTSINGLE\tCLUMPS\tOPTDUPL\tDEGEN\tQUALFILTERED\tTOOSHORT\tEMPTY" > $OUTPUTFOLDER/trimming.stats

for (( i = 1 ; i < $N+1 ; i++))
 do
SAMPLE=$(cat $LISTtoPROCESS | sed -n $i'p')
#collect stats

#READS=$(cat $OUTPUTFOLDER/$SAMPLE.optical_duplicates.log | grep  "Reads In" | cut -d":" -f2 | sed "s/ //g")
#READSIN=$(printf "%0.0f" $(bc -l <<< scale=6\;($READS / 2)) | sed "s/,/./g")

#136495552 reads processed; of these:

CLUMPS=1
OPTDUPL=1
#CLUMPS=$(cat $OUTPUTFOLDER/$SAMPLE.optical_duplicates.log | grep  "Clumps Formed" | cut -d":" -f2 | sed "s/ //g")
#OPTDUPL=$(cat $OUTPUTFOLDER/$SAMPLE.optical_duplicates.log | grep  "Duplicates Found" | cut -d":" -f2 | sed "s/ //g")

#SAMPLENAME=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "Input file" | rev | cut -d"/" -f1 | cut -d"." -f3-  | cut -d"_" -f2-| rev)

READPAIRSbefore=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "reads processed" | cut -d" " -f1)
READS=$READPAIRSbefore
READSIN=$READPAIRSbefore
#READSIN=$(printf "%0.0f" $(bc -l <<< scale=6\;($READS / 2)) | sed "s/,/./g")


DEGEN=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "degenerative reads" | cut -d"(" -f1 | sed "s/ //g")
QUALFILTERED=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "reads filtered out by quality control" | cut -d"(" -f1 | sed "s/ //g")
TOOSHORT=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "short reads filtered out after trimming by size control" | cut -d"(" -f1 | sed "s/ //g")
EMPTY=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "empty reads filtered out after trimming by size control" | cut -d"(" -f1 | sed "s/ //g")
READPAIRSafter=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "reads available; of these" | cut -d"(" -f1 | sed "s/ //g")
READPAIRSafterTRIMMED=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  " trimmed reads available after processing" | cut -d"(" -f1 | sed "s/ //g")
READPAIRSafterUNTRIMMED=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  " untrimmed reads available after processing" | cut -d"(" -f1 | sed "s/ //g") 
READPAIRSafterPC=$(cat $OUTPUTFOLDER/$SAMPLE.skewer.log | grep  "reads available; of these" | cut -d"(" -f2- | cut -d")" -f1| sed "s/ //g")


### commented out and given arbitrary value since sickle is not run here
SICKLEOUTPAIRED=1
SICKLEOUTSINGLE=1
SICKLEPC=100
#SICKLEOUTPAIRED=$(cat $OUTPUTFOLDER/strict/$SAMPLE.sickle.log | grep  "paired records kept" | cut -d"(" -f2 | cut -d" " -f1 | sed "s/ //g")
#SICKLEOUTSINGLE=$(cat $OUTPUTFOLDER/strict/$SAMPLE.sickle.log | grep  "single records kept" | cut -d":" -f2 | cut -d"(" -f1 | sed "s/ //g")

###SICKLEPC=$(echo $SICKLEOUTPAIRED * 100 / $READSIN | bc -l)
#SICKLEPC=$(printf "%0.2f" $(bc -l <<< scale=6\;($SICKLEOUTPAIRED * 100 / $READSIN)) | sed "s/,/./g")

echo "$SAMPLE ___ $READPAIRSafterPC __ ( $SICKLEPC% ) __ remain"
echo "$SAMPLE\t$READSIN\t$READPAIRSbefore\t$READPAIRSafter\t$READPAIRSafterPC\t$READPAIRSafterTRIMMED\t$READPAIRSafterUNTRIMMED\t$SICKLEOUTPAIRED\t$SICKLEPC\t$SICKLEOUTSINGLE\t$CLUMPS\t$OPTDUPL\t$DEGEN\t$QUALFILTERED\t$TOOSHORT\t$EMPTY" >> $OUTPUTFOLDER/trimming.stats

done

#plot stuff with R

### $PLOTSCRIPT $OUTPUTFOLDER $file
$PLOTSCRIPT $OUTPUTFOLDER trimming.stats

mv -f $OUTPUTFOLDER/*.pdf $OUTPUTFOLDER/QC/
mv -f $OUTPUTFOLDER/trimming.stats $OUTPUTFOLDER/QC/




################################################################################################
################################################################################################v
#### as single steps if wanted
#filter optical duplicates with clumpify, then filter with skewer, rename output and move to final folder
#cat $LISTtoPROCESS | parallel -j 1 "clumpify.sh in1=$INPUTFOLDER/{}'_1'.fastq.gz in2=$INPUTFOLDER/{}'_2'.fastq.gz out1=$TMPFOLDER/{}.temp.R1.fastq.gz out2=$TMPFOLDER/{}.temp.R2.fastq.gz dedupe optical dist=100 passes=1 subs=0 k=31 spantiles=f 2> >(tee $OUTPUTFOLDER/{}.optical_duplicates.log >&2)"

#filter optical duplicates with clumpify
#cat $LISTtoPROCESS | parallel -j 1 "clumpify.sh in1=$INPUTFOLDER/{}'_1'.fastq.gz in2=$INPUTFOLDER/{}'_2'.fastq.gz out1=$OUTPUTFOLDER/{}.temp.R1.fastq.gz out2=$OUTPUTFOLDER/{}.temp.R2.fastq.gz dedupe optical dist=100 passes=1 subs=0 k=31 spantiles=f 2> >(tee $OUTPUTFOLDER/{}.optical_duplicates.log >&2)"

#skewer, main filtering
#cat $LISTtoPROCESS | parallel -j 1 "skewer -m pe $INPUTFOLDER/{}.temp.R1.fastq.gz $INPUTFOLDER/{}.temp.R2.fastq.gz -x $FWD_ADAPTERS -y $REV_ADAPTERS --mean-quality $AVG_QUAL_THRESHOLD --end-quality $END_QUAL_THRESHOLD -l $L -n yes -r 0.1 -z --format auto -t $CPUs -o $OUTPUTFOLDER/{} && echo {} && mv $OUTPUTFOLDER/{}'-trimmed'.log $OUTPUTFOLDER/{}.skewer.log"
#cat $LISTtoPROCESS | head -n2 | parallel -j 1 "skewer -m pe $OUTPUTFOLDER/{}.temp.R1.fastq.gz $OUTPUTFOLDER/{}.temp.R2.fastq.gz -x $FWD_ADAPTERS -y $REV_ADAPTERS --mean-quality $AVG_QUAL_THRESHOLD --end-quality $END_QUAL_THRESHOLD -l $L -n yes -r 0.1 -z --format auto -t $CPUs -o $OUTPUTFOLDER/{} && echo {} && mv $OUTPUTFOLDER/{}'-trimmed'.log $OUTPUTFOLDER/{}.skewer.log"

# filter orphans with sickle ### skipped now (sickle seems to filter out way too much, orphans should have been filtered by skewer already
#SICKLE ; remove lowQ reads (and pair)#minL=50bp, minQ=15(sliding window), 4min per sample
#cat $LISTtoPROCESS | head -n2  | parallel -j $CPUs "sickle --mode pe -t sanger --gzip-output -q $AVG_QUAL_THRESHOLD -l $L -f $OUTPUTFOLDER/{}'-trimmed-pair1'.fastq.gz -r $OUTPUTFOLDER/{}'-trimmed-pair2'.fastq.gz -o $OUTPUTFOLDER/{}.R1.fastq -p $OUTPUTFOLDER/{}.R2.fastq -s $OUTPUTFOLDER/{}.singletons.fastq 1> >(tee $OUTPUTFOLDER/{}.sickle.log >&1)"

#renaming
#cat $LISTtoPROCESS | parallel -j 1 "mv $OUTPUTFOLDER/{}'-trimmed-pair1'.fastq.gz $OUTPUTFOLDER/{}.R1.fastq.gz && mv $OUTPUTFOLDER/{}'-trimmed-pair2'.fastq.gz $OUTPUTFOLDER/{}.R2.fastq.gz && mv $OUTPUTFOLDER/{}'-trimmed'.log $OUTPUTFOLDER/{}.skewer.log"

####
#fastqc
#mkdir $OUTPUTFOLDER/fastqcs
#ls -1 $OUTPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel -j 1 "fastqc --noextract -f fastq --threads 10 --nogroup -o $OUTPUTFOLDER/fastqcs/ $OUTPUTFOLDER/{}.fastq.gz && echo {}"

#multiqc
#mkdir $OUTPUTFOLDER/multiqc
#multiqc -o $OUTPUTFOLDER/multiqc/ $OUTPUTFOLDER/fastqcs/

#compress ##see alternativ with pigz (faster)
#tar -czf $OUTPUTFOLDER/fastqcs.tar.gz $OUTPUTFOLDER/fastqcs/





