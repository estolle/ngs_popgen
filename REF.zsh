#!/bin/zsh
############################################################################
# prep a new Ref.fa
# Eckart Stolle, Sept 2017
############################################################################
## ./REF.zsh INPUTFASTA 
#chmod -R 755 $HOME/scripts/REF.zsh

#REF="$HOME/ref/ref.fa.gz"


############################################################################

if [ $# -ne 1 ]; then
    echo $0: usage: ./REF.zsh INPUTFASTA
	echo "\nINPUTFASTA - new reference fasta (eg. .fa or .fas or .fna or .fasta) file (fix header names and uppercase/lowercase beforehand if needed"
    exit 1
fi

### check if progz are installed (numsum is optional, then comment/uncomment below (GZ/NN calcs)
which seqtk
which bedtools
which bwa
which samtools
which numsum
which python
which trf
which parallel

#mkdir -p $HOME/scripts/
CURRENTF=$PWD
if [[ ! -f $HOME/scripts/TRFdat_to_bed.py ]]; then
echo "DL script to convert TRFdat to bed"
cd $HOME/scripts
git clone https://github.com/hdashnow/TandemRepeatFinder_scripts/
mv $HOME/scripts/TandemRepeatFinder_scripts/TRFdat_to_bed.py $HOME/scripts/
mv $HOME/scripts/TandemRepeatFinder_scripts/TRFdat_to_txt.py $HOME/scripts/
chmod 777 $HOME/scripts/TRFdat_to_bed.py $HOME/scripts/TRFdat_to_txt.py
cd $CURRENTF
fi



#####################
#set/get variables
REF=$1

REFFOLDER=$(echo $REF | rev | cut -d"/" -f2- | rev)



if [[ ! -f $REF.amb ]]; then
echo "reference index not existing or not complete, creating index with bwa"
    bwa index $REF
fi


if [[ -f $(echo $REF | rev | cut -d"." -f2- | rev).gz ]]; then
echo "reference fasta gzipped"
  if [[ ! -f $(echo $REF | rev | cut -d"." -f2- | rev) ]]; then
    echo "unpacked ref file not found, unpacking"
    gunzip -k $REF
    DELETE1="$(echo $REF | rev | cut -d"." -f2- | rev)"
  fi
  echo "unpacked ref found"
  REF=$(echo $REF | rev | cut -d"." -f2- | rev)
  echo $REF
fi



##### name modifications if necessary
#zcat GCF_000188075.1_Si_gnG_genomic.fna.gz | fasta_formatter -w 0 > gng-mod1.fa
## sed -e 's/^\(>[^[:space:]]*\).*/\1/' $REF | fasta_formatter -w 0 > $REF.nodescriptors.fa
## cat $REF | fasta_formatter -w 0 > $REF.fixed

## convert lower case (repetitive/low complexity) to N and recreate Nbed
if [[ ! -f $REF.lowcmplx_to_N.fa ]]; then
echo "converting lower case into N (saved as new file)"
   cat $REF | sed -e '/^>/! s/[[:lower:]]/N/g' | seqtk seq -U > $REF.lowcmplx_to_N.fa
   COMPRESS1="$REF.lowcmplx_to_N.fa"
fi


## uppercase REF (-C drop fasta comments)
if [[ ! -f $REF.capitalized.fa ]]; then
echo "converting lower case into Uppercase (saved as new file)"
   cat $REF | seqtk seq -U -C > $REF.capitalized.fa
   COMPRESS2="$REF.capitalized.fa"
fi

#samtools index
if [[ ! -f $REF.fai ]]; then
echo "creating samtools index"
   samtools faidx $REF
fi

## create some simplified scaffold lists with length, also as regions-file (as input eg for samtools/SNPcall)
cat $REF.fai | cut -f1,2 > $REF.genometable
cat $REF.fai | cut -f1,2 > $REF.length
cat $REF.length | awk -F "\t" 'BEGIN {OFS = "\t";} {print $1,":0-",$2;}' | tr '\t' 'x' | sed "s/x//g" > $REF.region
cat $REF.length | awk -F "\t" 'BEGIN {OFS = "\t";} {print $1,"0",$2;}' > $REF.bed

## get some STATS
echo $REF.fai | rev | cut -d"/" -f1 | rev 1> >(tee > $REF.stats >&2)
echo "Number of Contigs and Scaffolds "$(cat $REF.fai | wc -l) 1> >(tee >> $REF.stats >&2)
echo "Total length of Contigs and Scaffolds "$(cat $REF.fai | numsum -c -x2) 1> >(tee >> $REF.stats >&2)
echo "STATS \n"$(cat $REF.fai | cut -f2 | sort -n | perl -ne 'chomp(); push(@contigs,$_);$total+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;if($sum>=$total*1.0){print "TOTAL: $total\n";exit;} ;}}' -) 1> >(tee >> $REF.stats >&2)
echo $(cat $REF.fai | cut -f2 | sort -n | perl -ne 'chomp(); push(@contigs,$_);$total+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;if($sum>=$total*0.5){print "N50 : $L\n";exit;} ;}}' -) 1> >(tee >> $REF.stats >&2)
echo $(cat $REF.fai | cut -f2 | sort -n | perl -ne 'chomp(); push(@contigs,$_);$total+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;if($sum>=$total*0.9){print "N90 : $L\n";exit;} ;}}' -) 1> >(tee >> $REF.stats >&2)

#total, n50, n90
#perl -e 'my ($len,$total)=(0,0);my @x;while(<>){if(/^[\>\@]/){if($len>0){$total+=$len;push@x,$len;};$len=0;}else{s/\s//g;$len+=length($_);}}if ($len>0){$total+=$len;push @x,$len;}@x=sort{$b<=>$a}@x; my ($count,$half)=(0,0);for (my $j=0;$j<@x;$j++){$count+=$x[$j];if(($count>=$total/2)&&($half==0)){print "TOTAL: $total\n";print "N50: $x[$j]\n";$half=$x[$j]}elsif($count>=$total*0.9){print "N90: $x[$j]\n";exit;}}' $REF 1> >(tee >> $REF.stats >&2)
## N50 oneliner (for a list of lengths)
#cat $REF.fai | cut -f2 | sort -n | perl -ne 'chomp(); push(@contigs,$_);$total+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;if($sum>=$total*0.5){print "TOTAL: $total\nN50 : $L\n";exit;} ;}}' -

## genomesize and Ns in the genome

#genomesize
## comment the numsum solution for GZ=... and uncomment the awk solution if you cannot install numsum
GZ=$(cat $REF.fai | numsum -c -x2) && echo $GZ
#GZ=$(cat $REF.fai | awk '{ sum+=$2} END {print sum}') && echo $GZ
GZ2=$(printf "%0.2f" $(bc -l <<< scale=6\;($GZ / 1000000))| sed "s/,/./g") && echo $GZ2

#number of N
## comment the numsum solution for GZ=... and uncomment the awk solution if you cannot install numsum
NN=$(seqtk comp $REF | cut -f9 | numsum -c -x1) && echo $NN
#NN=$(seqtk comp $REF | cut -f9 | awk '{ sum+=$1} END {print sum}') && echo $NN
NFRAC=$(printf "%0.2f" $(bc -l <<< scale=6\;(100-$NN*100/$GZ))| sed "s/,/./g") && echo $NFRAC
NFRAC2=$(printf "%0.2f" $(bc -l <<< scale=6\;($NN*100/$GZ))| sed "s/,/./g") && echo $NFRAC2

echo "Genomesize(MB): "$GZ2 1> >(tee >> $REF.stats >&2)
echo "Ns: "$NN 1> >(tee >> $REF.stats >&2)
echo "fraction of known bases (%): "$NFRAC 1> >(tee >> $REF.stats >&2)
echo "fraction of Ns (%): "$NFRAC2 1> >(tee >> $REF.stats >&2)


### make bedfile for Ns
python $HOME/scripts/fastaN2Bed.py $REF > $REF.NNN.bed
#python $HOME/scripts/fastaN2Bed.py $COMPRESS1 > $COMPRESS.NNN.bed

#NNNBED="$REF.NNN.bed"
#cat $REF.NNN.bed | cut -f 1 | sort | uniq > $REF.NNN.SCFs.lst
#mkdir /home/estolle/scratch/Soli/gnG/gnG_20161010-perscf-beds
#cat $REF.NNN.SCFs.lst | parallel -j $CPUs "grep {} $REF.NNN.bed > gnG_20161010-perscf-beds/{}.bed"
#cat /home/estolle/scratch/Soli/gnG/gnG_20161010.fa.NNN.bed | grep -v -E 'NC_014672|AF243435|AF243436|CP001391|AM999887' > $REF.NNN.bed

## find SSRs with TRF
################### SSRs step by step
echo "running tandem repeats finder (might take a while)"
echo "18\n36\n45\n48\n60" | parallel -j 5 --no-notice "echo 'running trf $REF 3 7 7 80 10 {} 5 -h -d'"
sleep 2s
echo "18\n36\n45\n48\n60" | parallel -j 5 --no-notice "trf $REF 3 7 7 80 10 {} 5 -h -d"

#trf $REF 3 7 7 80 10 18 5 -h -d
#trf $REF 3 7 7 80 10 36 5 -h -d
#trf $REF 3 7 7 80 10 45 5 -h -d
#trf $REF 3 7 7 80 10 48 5 -h -d
#trf $REF 3 7 7 80 10 60 5 -h -d

trf $REF 3 5 5 80 10 35 3 -h -d

mv *.3.7.7.80*.dat $REFFOLDER
mv *.3.5.5.80*.dat $REFFOLDER

DIR=$PWD
cd $REFFOLDER

TRFFILES=$(ls -1 *.3.7.7.80*.dat | rev | cut -d"." -f2- | rev)
N=$(echo $TRFFILES | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
SSRFILE=$(echo $TRFFILES | sed -n $i'p')
echo $SSRFILE
python $HOME/scripts/TRFdat_to_bed.py --dat $SSRFILE.dat --bed $SSRFILE.TRF.bed
gzip -9 $SSRFILE.dat
cat $SSRFILE.TRF.bed | sort -n | sortBed | awk -F "\t" 'BEGIN {OFS = "\t";} {$5=sprintf("%.0f",(($3-$2)/length($4)))} {print $1,$2,$3,"SSR_"$4"_"$5,(length($4));}' | awk -v "samplenumber=$i" '(($5) == samplenumber)' > $SSRFILE.SSRs.$i.nt.bed
done

SSRFILE=$(ls -1 *.3.5.5.80.10.35.3.dat | head -n1)
python $HOME/scripts/TRFdat_to_bed.py --dat $SSRFILE --bed $SSRFILE.TRF.bed
gzip -9 $SSRFILE
cat $SSRFILE.TRF.bed | sort -n | sortBed | awk -F "\t" 'BEGIN {OFS = "\t";} {$5=sprintf("%.0f",(($3-$2)/length($4)))} {print $1,$2,$3,"SSR_"$4"_"$5,(length($4));}' |awk '(($5) < 4)' > $SSRFILE.lowcomplx.bed

#fuse some SSRs bed (most relevant are 1-4nt, fornat/sort/hipster format
cat $REF.*.SSRs.1.nt.bed $REF.*.SSRs.2.nt.bed $REF.*.SSRs.3.nt.bed $REF.*.SSRs.4.nt.bed | awk -F "\t" 'BEGIN {OFS = "\t";} {print $1,$2-1,$3,$4,$5;}' | sort -n | sortBed > $REF.SSRs.bed
cat $REF.*.lowcomplx.bed | sort -n | sortBed > $REF.lowcmplx.bed
cat $REF.SSRs.bed | awk -F "\t" 'BEGIN {OFS = "\t";} {$5=sprintf("%.0f",(($3-$2)/length($4)))} {print $1,$2,$3,$5,$4;}' > $REF.SSRs.hipstr.bed

#create bedfile extended 5bp around each locus
bedtools slop -i $REF.SSRs.bed -g $REF.genometable -b 5 | sortBed > $REF.SSRs.ext5bp.bed
bedtools slop -i $REF.lowcmplx.bed -g $REF.genometable -b 5 | sortBed  > $REF.lowcmplx.ext5bp.bed

## tbd
# use bedtools to intersect N's and N's from the lower case (low complexity)

echo "done, cleaning up"

### compress & cleanup
#gzip -9 $COMPRESS1
#gzip -9 $COMPRESS2

#pigz -9 -p 20 $COMPRESS1
#pigz -9 -p 20 $COMPRESS2

pigz -9 -p 20 $REF.capitalized.fa
pigz -9 -p 20 $REF.lowcmplx_to_N.fa

#if [[ -f $DELETE1 ]]; then
#echo "deleting temporary ref file"
# rm $DELETE1
#fi

echo "reference fasta prepared"

