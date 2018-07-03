#!/bin/zsh

if [ $# -ne 1 ]; then
    echo $0: usage: $HOME/scripts/mosdepth.zsh INPUTBAM 
	echo "\nINPUTBAM: bam file"
    exit 1
fi

INPUT=$1
#INPUT="bam/Sol_U14-1-littleb-p_clade7B_712505-batch2.bam"
BAM=$(echo $INPUT | rev | cut -d"/" -f 1 | cut -d"." -f 2 | rev)
FOLDER=$(echo $INPUT | rev | cut -d"/" -f2- | rev)
mkdir -p $FOLDER/bamdepth

###
#mosdepth --threads 4 --by 5000 --mapq 55 $BAM.bam | bgzip -f -@4 -c /dev/stdin > $BAM.depth.500.bed.gz
#AVERAGE=$(zcat $BAM.depth.500.bed.gz | awk -v N=4 '{ sum += $N } END { if (NR > 0) print sum / NR }')
#mosdepth --threads 4 --by 10 --mapq 0 $BAM.bam | bgzip -f -@4 -c /dev/stdin > $BAM.depth.10LQ.bed.gz
#mosdepth --threads 4 --by 10 --mapq 55 $BAM.bam | bgzip -f -@4 -c /dev/stdin > $BAM.depth.10HQ.bed.gz

#mosdepth --threads 20 --by 10 --chrom KQ438548.1 --flag 4 test.Q01.lq10 bam/Q01.bam > test.Q01.lq10.bed
#cat test.Q01.lq10.bed | grep KQ438548

mosdepth --no-per-base --threads 90 --by 10 --mapq 0 --flag 4 $FOLDER/bamdepth/$BAM.depth.LQ $FOLDER/$BAM.bam
echo "$mosdepth 1 (incl LowQ) done"
mosdepth --no-per-base --threads 90 --by 10 --mapq 55 --flag 3844 $FOLDER/bamdepth/$BAM.depth.HQ $FOLDER/$BAM.bam
echo "$mosdepth 2 (HighQ) done"
mosdepth --no-per-base --threads 90 --by 5000 --mapq 55 --flag 3844 $FOLDER/bamdepth/$BAM.depth.HQ5k $FOLDER/$BAM.bam
echo "$mosdepth 3 (5kb HighQ) done"

#AVERAGE=$(zcat $FOLDER/bamdepth/$BAM.depth.HQ5k.regions.bed.gz | awk '$4 > 0 && $4 < 1000' | awk -v N=4 '{ sum += $N } END { if (NR > 0) print sum / NR }')
AVERAGE_STDEV=$(zcat $FOLDER/bamdepth/$BAM.depth.HQ5k.regions.bed.gz | awk '$4 > 0 && $4 < 1000' | awk -v N=4 '{FS=OFS="\t";}{x+=$N;y+=$N^2}END{print x/NR,sqrt(y/NR-(x/NR)^2)}')
AVERAGE=$(echo $AVERAGE_STDEV | cut -f1)
STDEV=$(echo $AVERAGE_STDEV | cut -f2)

##AVERAGE=$(grep "$BAM" $FOLDER/$BAM.depth.sums | cut -f2)
##AVERAGE=$(zcat $FOLDER/bamdepth/$BAM.depth.500.bed.gz | awk -v N=4 '{ sum += $N } END { if (NR > 0) print sum / NR }')
#AVERAGE=$(zcat $FOLDER/bamdepth/$BAM.depth.500.bed.gz | awk '$4 > 0 && $4 < 1000' | awk -v N=4 '{ sum += $N } END { if (NR > 0) print sum / NR }')

#zcat Q01.depth.LQ.regions.bed.gz | awk -v STDEV=$STDEV -v AVERAGE=$AVERAGE '$4 > (AVERAGE+2*STDEV) && $4 > 4' | bedtools merge -d 20 -c 4 -o mean,median,min,max -delim "\t" -i - | awk 'BEGIN{FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3-$2)}'

#high coverage = 2.0x the average coverage, but mind cov > 4 (meaning average cov of 50, min Cov 100)
zcat $FOLDER/bamdepth/$BAM.depth.LQ.regions.bed.gz | awk -v STDEV=$STDEV -v AVERAGE=$AVERAGE '$4 > (AVERAGE+2*STDEV) && $4 > 4 ' | bedtools merge -d 20 -c 4 -o mean,median,min,max -delim "\t" -i - | awk 'BEGIN{FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3-$2)}' | bgzip -f -@4 -c /dev/stdin > $FOLDER/bamdepth/$BAM.depth.LQ.regions.hicov.gz
ls $FOLDER/bamdepth/$BAM.depth.LQ.regions.hicov.gz
#zcat $FOLDER/bamdepth/$BAM.depth.LQ.per-base.bed.gz | awk -v AVERAGE=$AVERAGE '$4 > (AVERAGE*2.5) && $4 > 4 ' | bedtools merge -d 20 -c 4 -o mean,median,min,max -delim "\t" -i - | awk 'BEGIN{FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3-$2)}' | bgzip -f -@4 -c /dev/stdin > $FOLDER/bamdepth/$BAM.depth.LQ.per-base.hicov.gz
zcat $FOLDER/bamdepth/$BAM.depth.HQ.regions.bed.gz | awk '$4 == 1' | bedtools merge -d 1 -c 4 -o mean,median,min,max -delim "\t" -i - | awk 'BEGIN{FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3-$2)}' | bgzip -f -@4 -c /dev/stdin > $FOLDER/bamdepth/$BAM.depth.HQ.regions.locov.gz
ls $FOLDER/bamdepth/$BAM.depth.HQ.regions.locov.gz
zcat $FOLDER/bamdepth/$BAM.depth.HQ.regions.bed.gz | awk '$4 == 0' | bedtools merge -d 1 -c 4 -o mean,median,min,max -delim "\t" -i - | awk 'BEGIN{FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3-$2)}' | bgzip -f -@4 -c /dev/stdin > $FOLDER/bamdepth/$BAM.depth.HQ.regions.zerocov.gz
ls $FOLDER/bamdepth/$BAM.depth.HQ.regions.zerocov.gz

HICOV=$(zcat $FOLDER/bamdepth/$BAM.depth.LQ.regions.hicov.gz | awk '{sum += $8} END {print sum}')
LOCOV=$(zcat $FOLDER/bamdepth/$BAM.depth.HQ.regions.locov.gz | awk '{sum += $8} END {print sum}')
ZEROCOV=$(zcat $FOLDER/bamdepth/$BAM.depth.HQ.regions.zerocov.gz | awk '{sum += $8} END {print sum}')
echo "$BAM\t$AVERAGE\t$STDEV\t$HICOV\t$LOCOV\t$ZEROCOV"

echo "$BAM\t$AVERAGE\t$STDEV\t$HICOV\t$LOCOV\t$ZEROCOV" > $FOLDER/bamdepth/$BAM.depth.sums

#gzip $FOLDER/bamdepth/$BAM.depth.500.bed
#gzip $FOLDER/bamdepth/$BAM.depth.10LQ.bed
#rm -f $FOLDER/bamdepth/$BAM.depth.HQ.per-base.bed.gz $FOLDER/bamdepth/$BAM.depth.HQ.regions.bed.gz
#zcat $BAM.depth.10HQ.bed.gz | awk -v AVERAGE=$AVERAGE '$4 > (AVERAGE*2)' | bedtools merge -d 20 -c 4 -o mean,median,min,max -delim "\t" -i - | head -n100


