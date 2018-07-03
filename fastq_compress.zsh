### incomplete script

CPUs=40
ls /scratch/data/*/*.fastq | parallel -j1 "pigz -7 --processes $CPUs --stdout {} > {}.gz; echo {}"
ls /scratch/data/*/*.xz | rev | cut -d"." -f2- | rev | parallel -j1 "xz -d -T $CPUs --stdout {}.xz | pigz -7 --processes $CPUs --stdout > {}.gz; echo {}.xz"

ls /scratch/data/*/*.xz | rev | cut -d"." -f2- | rev | parallel -j1 "echo {}; rm -f {}.xz"
ls /scratch/data/*/*.fastq | parallel -j1 "echo {}; rm -f {}.fastq"


#########################

INPUTFOLDER="scatch/data/reads"
OUTPUTFOLDER="scratch/compressed_reads"
mkdir -p $OUTPUTFOLDER
#cp $HOME/progz/fqzcomp-4.6.tar.gz $OUTPUTFOLDER
cp /usr/local/src/archive/fqzcomp-4.6.tar.gz $OUTPUTFOLDER

DIRLIST=$(ls -d $INPUTFOLDER/* | rev | cut -d"/" -f 1 | rev)

N=$(echo $DIRLIST | wc -l | cut -d" " -f1)
echo "compressing $N folders for archiving"
for (( i = 1 ; i < $N+1 ; i++))
 do
FOLDER=$(echo $DIRLIST | sed -n $i'p')
echo $FOLDER
mkdir -p $OUTPUTFOLDER/$FOLDER
ls -1 $INPUTFOLDER/$FOLDER/*.fastq.gz | rev | cut -d"." -f 2- | cut -d"/" -f 1 | rev | parallel --no-notice -j 7 "echo {}; zcat $INPUTFOLDER/$FOLDER/{}.gz | fqz_comp -Q0 -s5+ -b -q3 > $OUTPUTFOLDER/$FOLDER/{}.fqz"
size $INPUTFOLDER
size $OUTPUTFOLDER
done




#INPUTFOLDER="/scratch/data/lane_349843"
#OUTPUTFOLDER="/scratch/data/compressed/lane_349843"
#mkdir -p $OUTPUTFOLDER
#
#ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 3- | cut -d"/" -f 1 | rev | parallel --no-notice -j 7 "echo {}; zcat $INPUTFOLDER/{}.fastq.gz | fqz_comp -Q0 -s5+ -b -q3 > #$OUTPUTFOLDER/{}.fastq.fqz"
#size $INPUTFOLDER
#size $OUTPUTFOLDER


## compress w xz
#ls -1 *tar.gz | rev | cut -d"." -f 3- | cut -d"/" -f 1 | rev | tail -n +4 | parallel --no-notice -j 1 "echo {}; zcat {}.tar.gz | xz -z -T 20 -9 --keep --stdout - > {}.tar.xz; echo done"



#compare  compression levels

#pigz -9 --keep --processes 30 --stdout $INPUTFILE1 > $INPUTFILE1.pigz9.gz
#pigz -9 --keep --processes 30 --stdout $INPUTFILE2 > $INPUTFILE2.pigz9.gz
#bgzip -@ 30 --stdout $INPUTFILE1 > $INPUTFILE1.bgzip.gz
#bgzip -@ 30 --stdout $INPUTFILE2 > $INPUTFILE2.bgzip.gz

###xz not parallel, takes ages, not suitable really, except if run parallel many times
#xz -z -9 --keep --stdout $INPUTFILE1 > $INPUTFILE1.xz9.gz
#xz -z -9 --keep --stdout $INPUTFILE2 > $INPUTFILE2.xz9.gz

#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.gz out2=$INPUTFILE2.clump.gz zl=9 pigz reorder
#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.bz2 out2=$INPUTFILE2.clump.bz2 zl=9 bzip2 reorder
#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.fqz out2=$INPUTFILE2.clump.fqz reorder

#time clumpify.sh in=$INPUTFILE1 out=$INPUTFILE1.clumpsingle.gz zl=9 pigz reorder
#time clumpify.sh in=$INPUTFILE2 out=$INPUTFILE2.clumpsingle.gz zl=9 pigz reorder

#time fqz_comp -Q0 -s5+ -b -q3 $INPUTFILE1 > $INPUTFILE1.fqz
#time fqz_comp -Q0 -s5+ -b -q3 $INPUTFILE2 > $INPUTFILE2.fqz




################
#pigz -9 --keep --processes 30 --stdout $INPUTFILE1 > $INPUTFILE1.pigz9.gz  1166,13s user 2,18s system 2736% cpu 42,697 total

#bgzip -@ 30 --stdout $INPUTFILE1 > $INPUTFILE1.bgzip.gz  302,08s user 6,75s system 2677% cpu 11,534 total
#bgzip -@ 30 --stdout $INPUTFILE2 > $INPUTFILE2.bgzip.gz  341,68s user 6,34s system 2651% cpu 13,127 total
#xz -z -9 --keep --stdout $INPUTFILE1 > $INPUTFILE1.xz9.gz  3467,56s user 3,02s system 99% cpu 57:51,42 total
#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.gz zl=9 pigz 2153,16s user 62,68s system 2413% cpu 1:31,83 total
#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.bz2 out2=$INPUTFILE2.clump.bz2 zl=9 bzip2 reorder Total time: 	346.191 seconds.
#clumpify.sh in1=$INPUTFILE1 in2=$INPUTFILE2 out1=$INPUTFILE1.clump.fqz out2=$INPUTFILE2.clump.fqz reorder  Total time: 	109.989 seconds.


#clumpify.sh in=$INPUTFILE1 out=$INPUTFILE1.clumpsingle.gz zl=9 pigz reorder  1232,32s user 50,43s system 1138% cpu 1:52,63 total
#Total time: 	111.351 seconds.

#clumpify.sh in=$INPUTFILE2 out=$INPUTFILE2.clumpsingle.gz zl=9 pigz reorder  1178,45s user 47,53s system 953% cpu 2:08,59 total
#Total time: 	127.223 seconds

#fqz_comp -Q0 -s5+ -b -q3 $INPUTFILE1 > $INPUTFILE1.fqz  166,24s user 1,86s system 146% cpu 1:54,91 total
#fqz_comp -Q0 -s5+ -b -q3 $INPUTFILE2 > $INPUTFILE2.fqz  171,64s user 2,01s system 155% cpu 1:51,88 total



### comparison of compressors

#xz9.gz not parallel, but maybe best because it doesnt change the files
#clump.bz2 346.191 seconds
#clumpify/fqz 109.989 seconds good compression
#.clump.fqz in paired order, fast
#fqz: good, but maybe not lossless, long even with multithread #in paired order, Ns quality (!) is set to #

#-rwxrwxr-x 1 es es 2,5G Sep 14 12:44 WTCHG_336855_005_1.fastq
#-rwxrwxr-x 1 es es 686M Sep 14 12:44 WTCHG_336855_005_1.fastq.gz
#-rw-rw-r-- 1 es es 651M Sep 14 16:53 WTCHG_336855_005_1.fastq.bgzip.gz
#-rw-rw-r-- 1 es es 615M Sep 14 16:51 WTCHG_336855_005_1.fastq.pigz9.gz

#-rw-rw-r-- 1 es es 524M Sep 14 21:36 WTCHG_336855_005_1.fastq.clumpsingle.gz
#-rw-rw-r-- 1 es es 508M Sep 14 21:19 WTCHG_336855_005_1.fastq.clump.gz

#-rw-rw-r-- 1 es es 486M Sep 14 17:51 WTCHG_336855_005_1.fastq.xz9.gz
#-rw-rw-r-- 1 es es 454M Sep 14 21:30 WTCHG_336855_005_1.fastq.clump.bz2
#-rw-rw-r-- 1 es es 400M Sep 14 21:34 WTCHG_336855_005_1.fastq.clump.fqz

#-rw-rw-r-- 1 es es 385M Sep 14 22:13 WTCHG_336855_005_1.fastq.fqz


###########################################
#-rwxrwxr-x 1 es es 2,5G Sep 14 13:03 WTCHG_336855_005_2.fastq
#-rwxrwxr-x 1 es es 808M Sep 14 13:03 WTCHG_336855_005_2.fastq.gz
#-rw-rw-r-- 1 es es 775M Sep 14 16:53 WTCHG_336855_005_2.fastq.bgzip.gz
#-rw-rw-r-- 1 es es 737M Sep 14 16:52 WTCHG_336855_005_2.fastq.pigz9.gz

#-rw-rw-r-- 1 es es 682M Sep 14 21:38 WTCHG_336855_005_2.fastq.clumpsingle.gz
#-rw-rw-r-- 1 es es 662M Sep 14 21:19 WTCHG_336855_005_2.fastq.clump.gz

#-rw-rw-r-- 1 es es 599M Sep 14 19:10 WTCHG_336855_005_2.fastq.xz9.gz
#-rw-rw-r-- 1 es es 580M Sep 14 21:30 WTCHG_336855_005_2.fastq.clump.bz2
#-rw-rw-r-- 1 es es 501M Sep 14 21:34 WTCHG_336855_005_2.fastq.clump.fqz

#-rw-rw-r-- 1 es es 486M Sep 14 22:15 WTCHG_336855_005_2.fastq.fqz


#2,5G .fastq
#808M .fastq.gz
#775M .fastq.bgzip.gz
#737M .fastq.pigz9.gz (gz compression in parallel with pigz, level 9)
#682M .fastq.clumpsingle.gz (clumpify on single read files)
#662M .fastq.clump.gz (clumpify on paired of reads)
#599M .fastq.xz9.gz (xz level 9)
#580M .fastq.clump.bz2 (clumpify on pairs and bz2)
#501M .fastq.clump.fqz (clumpify on pairs and fqzcomp)
#486M .fastq.fqz (fqzcomp)













