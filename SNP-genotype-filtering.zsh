#!/bin/zsh
############################################################################
# filter SNPs
# Eckart Stolle, Jan.2018
############################################################################
## script: SNP-genotype-filtering.zsh
## scripts processes 1 file, outputs filtered VCF
## use parallel to run it simultaneous on several files
## uses 2-5 compression threads
## makes a temporary folder in the RAM, thus watch the RAM usage

if [ $# -ne 3 ]; then
    echo $0: usage: ./SNP-genotype-filtering.zsh INPUT.raw.vcf.gz SAMPLEname OUTPUTFOLDER 
	echo "\nINPUT.raw.vcf.gz: vcf file to filter (bgzip-ed/tabix-indexed"
	echo "\nSAMPLEname: Sample Name, will be used for output VCF file"
	echo "\OUTPUTFOLDER: output folder for filtered VCF file"
    exit 1
fi

#set/get variables
VCFINPUT=$1
SAMPLE=$2
FILTEROUTFOLDER=$3
CPUs=10


# !! Check your VCFs genotype declaration
echo "genotype field format/abbreviation used here: GT:GQ:DP:AD:RO:QR:AO:QA:GL"



# make a temporary folder in the RAM
TF="/dev/shm/$SAMPLE"
mkdir -p $TF

# we are changing DIR into that folder to prevent bgzip mangeling all the temporary files created in parallel from different vcf's. those temp names are all the same, causing problems
ORIGFOLDER=$PWD
cd $TF

#### set limits for filtering
COVINCL=2
FRACTION_AO_RO=0.25

####################################################################################
#### coverage at which a SNP is not considered


#COVINCL=4

echo "reads coverage threshold: "$COVINCL



####################################################################################
# optional:
#### coverage at which a SNP is not considered
##
## automatically IncludedCoverage-1
##
##
COVEXCL=$(echo $COVINCL-1 | bc -l)
##
#echo "--> reads coverage excluded: "$COVEXCL
##


####################################################################################
#### lower limit for Fraction AO/RO (reads w alternative allele / reads w ref-allele)
##
#### eg 0.25 --> limit is 4 reads w alternativ allele vs. 16 reads w ref allele
#### for 1/3 put 0.33333333333333333333333333
##
## this fraction limit is inclusive
##   --> in case there are 4 reads ALT and 16 reads REF (4/16=0.25),
##       the Locus will be considered Heterozygous
##   --> in case there are 3 reads ALT and 15 reads REF (3/15=0.2),
##       the Locus is considered homozygous
##
#### such a cutoff is very helpful in the case of haploid/diploid samples to avoid spurious heterozygous calls
#### expected ALT/REF in a diploid locus ideally is 1:1, but there is typically some variability
#### highly PCRed samples can have igher variance
#### for haploids, a higher threshold might be beneficial (eg. 0.33333333)


#FRACTION_AO_RO=0.25
echo "threshold ALT/REF reads for SNP from which considered heterozygous: "$FRACTION_AO_RO



######################################################################################
##
# set other threshold automatically
##
FRACTION_AO_ROplus=$(echo $FRACTION_AO_RO+0.00000000000000000000000000000000001 | bc -l)
FRACTION_AO_ROminus=$(echo $FRACTION_AO_RO-0.00000000000000000000000000000000001 | bc -l)
##
echo  "--> faction ALT/REF considered heterzygous: "$FRACTION_AO_ROplus
echo  "--> faction ALT/REF considered homozygous: "$FRACTION_AO_ROminus
##

####################################################################################
# optional:
#### lower limit for Fraction AO/all (reads w alternative allele / all reads
##
#### optionally can be set, otherwise determined automatically
##
#### eg 0.2 --> limit is 4 reads w alternativ allele in 16 reads total --> 4readsALT vs 12readsREForOTHER=0.2
#### for 1/3 put 0.33333333333333333333333333
##
## this fraction limit is inclusive
##   --> in case there are 4 reads ALT and 16 reads in total (=12 reads not ALT) (ALT/REF: 4/12=0.20),
##       the Locus will be considered Heterozygous
##   --> in case there are 3 reads ALT and 15 reads in total (=12 reads not ALT) (ALT/REF: 3/12=0.25),
##       the Locus is considered Heterozygous
##   --> in case there are 3 reads ALT and 19 reads in total (=16 reads not ALT) (ALT/REF: 3/16=0.1875),
##       the Locus is considered Homozygous
##
#### such a cutoff is very helpful in the case of haploid/diploid samples to avoid spurious heterozygous calls
#### expected ALT/REF in a diploid locus ideally is 1:1, but there is typically some variability
#### highly PCRed samples can have igher variance
#### for haploids, a higher threshold might be beneficial (eg. 0.25)
#### this total threshold is somewhat redundant to the above threshold for the fraction of ALT/REF
####     but in the case of many alternative alleles or other reads at the locus, this might be helptful
##
##
FRACTION_total=$(echo "1/(1/$FRACTION_AO_RO+1)" | bc -l)
FRACTION_totalplus=$(echo $FRACTION_total+0.00000000000000000000000000000000001 | bc -l)
FRACTION_totalminus=$(echo $FRACTION_total-0.00000000000000000000000000000000001 | bc -l)
##
echo "threshold ALT/TotalReads for SNP from which considered heterozygous: "$FRACTION_total
echo "--> considered heterozygous: "$FRACTION_totalplus
echo "--> considered homozygous: "$FRACTION_totalminus
##
##

####################################################################################
# optional:
#### upper limit for Fraction AO/all (reads w alternative allele / all reads
##
## optionally set the threshold, otherwise will be det automatically from Fraction AO/all above
## (1-that fraction)
##
#### eg 0.75 --> limit is 12 reads w alternativ allele vs. 16 reads in total
##   --> 75% of reads are ALT allele
##   --> considered homozygous for ALT allele
##
## this fraction limit is inclusive, threshold == heterozygous, above == homozygous ALT
##   --> in case there are 12 reads ALT and 16 reads total (=4 reads ref) (12/16=0.75, 12ALT/4REF=3),
##       the Locus will be considered heterozygous, analogous for the lower limit
##   --> in case there are 12 reads ALT and 15 reads total (=3 reads ref) (12/15=0.8, 12ALT/5REF=4),
##       the Locus will be considered homozygous for ALT allele


FRACTIONuLIMIT=$(echo "1-$FRACTION_total" | bc -l)
FRACTIONuLIMITPLUS=$(echo $FRACTIONuLIMIT+0.00000000000000000000000000000000001 | bc -l)
FRACTIONuLIMITMINUS=$(echo $FRACTIONuLIMIT-0.00000000000000000000000000000000001 | bc -l)
echo "threshold ALT/all_reads for SNP from which considered homozygous for ALT allele: "$FRACTIONuLIMIT


####################################################################################
# thresholds summary:
echo
echo "thresholds summary: \nCoverage considered: "$COVINCL" X"
echo "homozygous (REF allele): <"$FRACTION_total
echo "homozygous (ALT allele): >"$FRACTIONuLIMIT
echo "heterozygous: "$FRACTION_total" to "$FRACTIONuLIMIT



####################################################################################

#### if run from different folder (to run in parallel but avoid bgzip mixup, as it is placing the same named tmp files in current dir)

#store threshold in the OUTPUTFOLDER
echo $FOLDERNAME > $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo $OUTPUTFOLDER >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo $FILTEROUTFOLDER >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "COVINCL\t"$COVINCL >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "COVEXCL\t"$COVEXCL >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_AO_RO\t"$FRACTION_AO_RO >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_AO_ROplus\t"$FRACTION_AO_ROplus >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_AO_ROminus\t"$FRACTION_AO_ROminus >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_total\t"$FRACTION_total >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_totalplus\t"$FRACTION_totalplus >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTION_totalminus\t"$FRACTION_totalminus >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTIONuLIMIT\t"$FRACTIONuLIMIT >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTIONuLIMITPLUS\t"$FRACTIONuLIMITPLUS >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt
echo "FRACTIONuLIMITMINUS\t"$FRACTIONuLIMITMINUS >> $FILTEROUTFOLDER/filter.vcf.thresholds.txt


##### get header
echo "$SAMPLE"
tabix -H $VCFINPUT | sed "s#ID=DPR,#ID=AD,#g" > $TF/$SAMPLE.filtered.header
echo "header extracted"

#change GT field DPR to AD for compatibility
##FORMAT=<ID=DPR,Number=A,Type=Integer,Description="Number of observation for each allele">
#sed "s#ID=DPR,#ID=AD,#g"

# Inputlines
LINESBEFORE=$(zcat $VCFINPUT | grep -v "#" | wc -l)
#2740291 with header
#2734578 without header

#cat $SAMPLE.vcf | vcffilter -g " DP = 2 " | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$"
## example of a lowcov
#Group10.10	0|1:20.8169:2:1,1:1:33:1:33:-2.69707,0,-2.69707
 #cat $SAMPLE.vcf | vcffilter -g " DP < $COVINCL " | grep -P Group10.10"\t"332420
#0|1:20.8169:2:1,1:1:33:1:33:-2.69707,0,-2.69707
 #cat $SAMPLE.vcf | vcffilter -g " DP > $COVEXCL " | grep -P Group10.10"\t"332420
#. #needs to be filtered
 #cat $SAMPLE.vcf | vcffilter -sg " DP > $COVEXCL " | grep -P Group10.10"\t"332420
# no GT:GQ:DP:AD:RO:QR:AO:QA:GL
 #cat $SAMPLE.vcf | vcffilter -g " DP > $COVEXCL " | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" | grep -P Group10.10"\t"332420

# missing data examples
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.:.:.:.:.:.:.:.:.
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.|.:8.54444e-12:1:0,0:0:0:0:0:0,0,0
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.|.:6.044e-08:2:0,0,0:0:0:0,0:0,0:0,0,0,0,0,0
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.|0:3.25922e-06:12:5,0,3:5:175:0,3:0,111:-7.54918,-9.05433,-23.2934,0,-14.2496,-13.3465
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.|1:0.00666328:7:0,5:0:0:5:157:-14.0816,-1.50515,0
#GT:GQ:DP:AD:RO:QR:AO:QA:GL	.

#works
 # vcffilter -g " DP > $COVEXCL " $VCFINPUT | tee >( grep -P 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\.\|\.|\.\/\.|\.:\.:\.:\.:\.' > test5) > $SAMPLE.vcf
#2740292
#221846
#cat $SAMPLE.vcf | grep -c -P 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\.\|\.|\.\/\.|\.:\.:\.:\.:\.'
#221846
#cat $SAMPLE.vcf | grep -cP 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$'
#221836

#cat test5 | sed "s#.:.:.:.:.:.:.:.:.#.|.#g" | sed "s#\t\.\$#\t.|.#g" | sed "s#\.\/\.#.|.#g" | cut -f10 | sort | uniq
#.|.
#.|.:0.000427112:4:0,0:0:0:0:0:0,-0.647376,-0.0453156
#.|.:1.68251e-07:4:0,0,0:0:0:0,0:0,0:0,0,0,0,0,0
#.|.:3.68379e-09:4:0,0,0,0:0:0:0,0,0:0,0,0:0,-0.647618,-0.0455579,0,-0.647618,0,0,-0.647618,0,0
#.|.:4.36251e-09:4:0,0,0:0:0:0,0:0,0:0,0,0,0,0,0
#.|.:5.88698e-09:4:0,0,0,0:0:0:0,0,0:0,0,0:0,0,0,0,0,0,0,0,0,0

###### 1. take out missing data

#! optional to build in (done above already) cat $SAMPLE.vcf | sed -r "s#GT:GQ:DP:AD:RO:QR:AO:QA:GL\t(.)\/(.)#GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\1|\2#g"
# to set everything to phased, later only the true heterozygous data can be set to unphased
#sed "s#\.\/\.#.|.#g"

##########################################   DATA-input, filter by coverage/missing
######  take out missing & low cov data
##########################################
# extended pipe with tee
vcffilter -g " DP > $COVEXCL " $VCFINPUT | sed "s#GT:GQ:DP:DPR:RO:QR:AO:QA:GL#GT:GQ:DP:AD:RO:QR:AO:QA:GL#g" | tee >( grep -P 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\t\.\|\.|\t\.\/\.|\.:\.:\.:\.:\.' |\
sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.\$#\t.|.#g" | sed -r "s#\t\.\|\.\$#\t.|.#g" | sed -r "s#\t\.\|\.#\t.#g" > $TF/$SAMPLE.filtered.missing) |\
grep -vP 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\t\.\|\.|\t\.\/\.|\.:\.:\.:\.:\.' > $TF/$SAMPLE.vcf


#cat $TF/$SAMPLE.vcf | sed "s#GT:GQ:DP:DPR:RO:QR:AO:QA:GL#GT:GQ:DP:AD:RO:QR:AO:QA:GL#g" | grep -P 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$'

NONMISSINGGENOTYPES=$(cat $TF/$SAMPLE.vcf | grep -v "#" | wc -l) && echo "nonmissing: "$NONMISSINGGENOTYPES
MISSINGGENOTYPES=$(cat $TF/$SAMPLE.filtered.missing | grep -v "#" | wc -l) && echo "missing: "$MISSINGGENOTYPES
FIRSTSUM=$(echo "$NONMISSINGGENOTYPES+$MISSINGGENOTYPES" | bc -l)
echo "before/after first filtering (sum): $LINESBEFORE $FIRSTSUM"
#2512732
#221846
#2734578

# if wanted, set all genotypes to phased, all heterozygous genotpes later need to be changed to "/" unphased
#cat $SAMPLE.vcf | grep -v "#" | cut -f10 | cut -d":" -f1 | sort | uniq
#.|0
#0|0
#0|1
#.|1
#1|0
#1|1

##### take out heterozygous/homozygous
#zgrep "0/1" ${VCFFILE}.$SAMPLE.vcf.gz | grep -v "#" | grep "NUMALT=2" | awk '(($5) ~ /,/)'

## Homzygous are clear
## Heterozygous split by having NUMALT=1 (biallelic) or not with NUMALT=1 (both with ALT/REF genotypes)
## the latter occurs with triallelic
## Heterozygous with ALT1/ALT2 genotyes are separated

rm -f $TF/$SAMPLE.filtered.homozygous
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz2
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz3
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.partiallymissing

cat $TF/$SAMPLE.vcf | grep -v "^#" | tee >(grep -P '0\|0|1\|1|2\|2|3\|3|0\/0|1\/1|2\/2|3\/3' >> $TF/$SAMPLE.filtered.homozygous) \
>(grep -v "NUMALT=1" | grep -P '1\|0|0\|1|2\|0|0\|2|0\|3|3\|0|0/1|1/0|0/2|2/0|0/3|3/0' >> $TF/$SAMPLE.tmp.htz2) \
>(grep -P '\.\|1|\.\|0|0\|\.|1\|\.|\.\|3|\.\|2|2\|\.|3\|\.' >> $TF/$SAMPLE.tmp.partiallymissing) \
>(grep -P '1\|2|2\|1|1\|3|3\|1|2\|3|3\|2|1/2|2/1|1/3|3/1|2/3|3/2' >> $TF/$SAMPLE.tmp.htz3) |\
grep "NUMALT=1" | grep -P '1\|0|0\|1|2\|0|0\|2|0\|3|3\|0|0/1|1/0|0/2|2/0|0/3|3/0' >> $TF/$SAMPLE.tmp.htz

rm -f $TF/$SAMPLE.vcf

HMZ=$(cat $TF/$SAMPLE.filtered.homozygous | grep -v "#" | wc -l) && echo "homozygous: "$HMZ
HTZ=$(cat $TF/$SAMPLE.tmp.htz | grep -v "#" | wc -l) && echo "biallelic heterozygous REF/ALT: "$HTZ
HTZ2=$(cat $TF/$SAMPLE.tmp.htz2 | grep -v "#"| wc -l) && echo "multiallelic heterozygous REF/ALT: "$HTZ2
HTZ3=$(cat $TF/$SAMPLE.tmp.htz3 | grep -v "#" | wc -l) && echo "multiallelic heterozygous ALT/ALT: "$HTZ3
PARTMISSING=$(cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#"  | wc -l) && echo "partially missing genotype (eg .|1): "$PARTMISSING


#homozygous: 1884185
#biallelic heterozygous REF/ALT: 627052
#multiallelic heterozygous REF/ALT: 0
#multiallelic heterozygous ALT/ALT: 0
#partially missing genotype (eg .|1): 1495

echo "homozygous loci extracted: "$HMZ
echo "heterozygous loci extracted: "$(echo "$HTZ+$HTZ2+$HTZ3" | bc -l)
#homozygous loci extracted: 1884185
#heterozygous loci extracted: 627052
SECONDSUM=$(echo "$HMZ+$HTZ+$HTZ2+$HTZ3+$PARTMISSING" | bc -l)
echo "homozygous+heterozygous: "$SECONDSUM" from: "$NONMISSINGGENOTYPES



#######################################   HTZ
#### analyze/filter the heterozygots
#######################################


##      filter false heterozygots, separate true heterozygots (can be set to missing missing value)
## NUMALT=1 and NUMALT=2/3 come in separate files, NUMALT=2/3: use only RO measurements (to avoid problem with the two ALT alleles)

# HTZ actually ALT homoz
#cat $TF/$SAMPLE.tmp.htz | vcffilter -g "RO / AO < $FRACTION_AO_ROplus" |\
#grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
#sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' >> $TF/$SAMPLE.filtered.htz.alt-hmz

# HTZ actually REF homoz
#cat $TF/$SAMPLE.tmp.htz | vcffilter -g "AO / RO < $FRACTION_AO_ROplus" |\
#grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
#sed -E -e 's#0\/1|1\/0|0\|1|1\|0|0\/2|2\/0|0\|2|2\|0#0|0#g' >> $TF/$SAMPLE.filtered.htz.ref-hmz

# True HTZ
#cat $TF/$SAMPLE.tmp.htz | vcffilter -g "RO / AO > $FRACTION_AO_ROplus" | vcffilter -g "AO / RO > $FRACTION_AO_ROplus" |\
#grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" >> $TF/$SAMPLE.filtered.htz.heterozygous

####### modify the AO field for calls which originally contained several alleles (AO field has several entries, causing vcffilter to make errors when using that field)
####### this stems from haplotypes found (and their alleles) even tho they were decomposed (and filtered for biallelic) (the genotype calls are correct, but the handling here for filtering not)
# works
#cat numbers | tr ',' '\n' | awk 'NR==1{max=$1} {if($1>max){max=$1}} END { print max }'
#does not work in parallel
#cat numbers | parallel echo {} | tr ',' '\n' | awk 'NR==1{max=$1} {if($1>max){max=$1}} END { print max }'
#cat numbers | parallel --pipe echo {} | tr ',' '\n' | awk \'NR==1{max=\$1} {if\(\$1>max\){max=\$1}} END { print max }\'
#cat SRR5115489.tmp.htz-b | parallel -k -j110 "echo {} | tr ',' '\n' | awk \'NR==1{max=\$1} {if ( \$1>max ){max=\$1 }} END { print max }\' " >> SRR5115489.tmp.htz-b-highest.test"

# same speed
#cat SRR5115489.tmp.htz-b | parallel -k -j110 "echo {} | tr ',' '\n' | numbound >> testnumbound" 
#cat SRR5115489.tmp.htz-b | parallel -k -j110 "echo {} | tr ',' '\n' | sort -rn | head -n1 >> testsort"
#seem lil faster, but its the sum
#cat SRR5115489.tmp.htz-b | parallel -k -j110 "echo {} | tr ',' '+' | bc -l >> testbc_sum"

#### working, but slow, see below using R and intermediate tables/lists
##  split vcf
#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f1-14 > $TF/$SAMPLE.tmp.htz-a
#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f15 > $TF/$SAMPLE.tmp.htz-b
#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f16- > $TF/$SAMPLE.tmp.htz-c
##  for the AO field, take the highest AO count, parallel calc the sum of AO counts (all this is relevant for several ALT alleles only, frequently left over after decomposition)
### runs v slow, no alternatives yet found
#cat $TF/$SAMPLE.tmp.htz-b | parallel -k -j10 "echo {} | sed 's#,#\n#g' | tee >(sort -rn | head -n1 >> $TF/$SAMPLE.tmp.htz-b-highest) | numsum -x1 >> $TF/$SAMPLE.tmp.htz-b-sum"
##  paste together again
##paste $TF/$SAMPLE.tmp.htz-b $TF/$SAMPLE.tmp.htz-b-highest | head -n10000
#cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz-modified.highest
#cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz-modified.sum
#paste -d":" $TF/$SAMPLE.tmp.htz-a $TF/$SAMPLE.tmp.htz-b-highest $TF/$SAMPLE.tmp.htz-c >> $TF/$SAMPLE.tmp.htz-modified.highest
#paste -d":" $TF/$SAMPLE.tmp.htz-a $TF/$SAMPLE.tmp.htz-b-sum $TF/$SAMPLE.tmp.htz-c >> $TF/$SAMPLE.tmp.htz-modified.sum
#rm -f $TF/$SAMPLE.tmp.htz-a $TF/$SAMPLE.tmp.htz-b $TF/$SAMPLE.tmp.htz-c $TF/$SAMPLE.tmp.htz-b-highest $TF/$SAMPLE.tmp.htz-b-sum

##############################################################

# working version, complete (uncommented) below

#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f1-14 > $TF/$SAMPLE.tmp.htz-a
#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f15 > $TF/$SAMPLE.tmp.htz-b
#cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f16- > $TF/$SAMPLE.tmp.htz-c
#cat $TF/$SAMPLE.tmp.htz-b | sed "s#,0##g" > $TF/$SAMPLE.fixgenotypes
#cat $TF/$SAMPLE.tmp.htz-b | sed "s#,0##g"  | awk -F ',' ' { print NF-1 } ' > $TF/$SAMPLE.commanumbers
#MAXCOMMA=$(cat commanumbers | numbound) && echo $MAXCOMMA
#cat $TF/$SAMPLE.commanumbers | awk -v MAXCOMMA="$MAXCOMMA" ' { print MAXCOMMA-$1}' > $TF/$SAMPLE.commas2add
#cat $TF/$SAMPLE.commas2add | sed "s#1#,#g" | sed "s#2#,,#g" | sed "s#3#,,,#g" | sed "s#4#,,,,#g" | sed "s#5#,,,,,#g" | sed "s#6#,,,,,,#g" | sed "s#7#,,,,,,,#g" | sed "s#8#,,,,,,,,#g" > $TF/$SAMPLE.commas2substituted
#paste -d "" $TF/$SAMPLE.fixgenotypes $TF/$SAMPLE.commas2substituted > $TF/$SAMPLE.fixgenotypes2
#-d'\0'
#-d ""

####### too slow version to process the numbers
#head -c 3 < /dev/zero | tr '\0' ','
#works, slow
#cat commas2add | parallel -k -j1 "head -c {} < /dev/zero | tr '\0' ',' >> commas; echo >> commas"
#works, slow
#cat commas2add | parallel -k -j1 'printf $(for ((i=1; i<={}; i++));do printf "%s" ",";done;printf "\n")"\n" >> commas2'

#N=4
#perl -e "print ',' x $N;"
#printf "%0.s," $(seq 1 $N)

#COMMA=$(printf "%0.s," $(seq 1 $N))
#echo $COMMA >> file
#echo $(printf "%0.s," $(seq 1 $N))
###########

##### fastest way, using an R script (see below)

###################################################################################################################
## #######################       R script to det the max number of any ALT AO number
# $HOME/scripts/maxnumber.R

#args <- commandArgs(trailingOnly = TRUE)
#print(args[1])
#inputfile <- (args[1])
###inputfile <- "/dev/shm/SRR5115527/fixgenotypes2"
#data <- read.csv(file=inputfile, header=FALSE, sep=",",colClasses="integer")
###max (data$V2,na.rm=TRUE)
#colMax <- apply(data, 1, function(x) max(x,na.rm=TRUE))
#data <- cbind(colMax, data)
#write.csv(data, file = paste0(inputfile,".csv") ,row.names=FALSE)
###################################################################################################################
###################################################################################################################

cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f1-14 > $TF/$SAMPLE.tmp.htz-a
cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f15 > $TF/$SAMPLE.tmp.htz-b
cat $TF/$SAMPLE.tmp.htz | grep -v "#" | cut -d":" -f16- > $TF/$SAMPLE.tmp.htz-c

cat $TF/$SAMPLE.tmp.htz-b | sed "s#,0##g" > $TF/$SAMPLE.tmp.htz.fixgenotypes
cat $TF/$SAMPLE.tmp.htz-b | sed "s#,0##g"  | awk -F ',' ' { print NF-1 } ' > $TF/$SAMPLE.tmp.htz.commanumbers
MAXCOMMA=$(cat $TF/$SAMPLE.tmp.htz.commanumbers | numbound) && echo $MAXCOMMA
cat $TF/$SAMPLE.tmp.htz.commanumbers | awk -v MAXCOMMA="$MAXCOMMA" ' { print MAXCOMMA-$1}' > $TF/$SAMPLE.tmp.htz.commas2add
cat $TF/$SAMPLE.tmp.htz.commas2add | sed "s#1#,#g" | sed "s#2#,,#g" | sed "s#3#,,,#g" | sed "s#4#,,,,#g" | sed "s#5#,,,,,#g" | sed "s#6#,,,,,,#g" | sed "s#7#,,,,,,,#g" | sed "s#8#,,,,,,,,#g" > $TF/$SAMPLE.tmp.htz.commas2substituted
paste -d "" $TF/$SAMPLE.tmp.htz.fixgenotypes $TF/$SAMPLE.tmp.htz.commas2substituted > $TF/$SAMPLE.tmp.htz.fixgenotypes2

Rscript $HOME/scripts/maxnumber.R $TF/$SAMPLE.tmp.htz.fixgenotypes2
cat $TF/$SAMPLE.tmp.htz.fixgenotypes2.csv | tail -n+2 | cut -d"," -f1 > $TF/$SAMPLE.tmp.htz.fixedcolumn

cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz-modified.highest
paste -d":" $TF/$SAMPLE.tmp.htz-a $TF/$SAMPLE.tmp.htz.fixedcolumn $TF/$SAMPLE.tmp.htz-c >> $TF/$SAMPLE.tmp.htz-modified.highest

rm -f $TF/$SAMPLE.tmp.htz.commanumbers $TF/$SAMPLE.tmp.htz.commas2add $TF/$SAMPLE.tmp.htz.commas2substituted $TF/$SAMPLE.tmp.htz.fixgenotypes2 $TF/$SAMPLE.tmp.htz.fixgenotypes
rm -f $TF/$SAMPLE.tmp.htz-a $TF/$SAMPLE.tmp.htz-b $TF/$SAMPLE.tmp.htz-c $TF/$SAMPLE.tmp.htz-b-highest $TF/$SAMPLE.tmp.htz.fixgenotypes2.csv


####

# Tee/Pipes for parallelization
rm -f $TF/$SAMPLE.filtered.htz.ref-hmz $TF/$SAMPLE.filtered.htz.alt-hmz $TF/$SAMPLE.filtered.htz.heterozygous
cat $TF/$SAMPLE.tmp.htz-modified.highest | tee \
>(vcffilter --or -g "RO / AO > $FRACTION_AO_RO | RO / AO = $FRACTION_AO_RO" |\
vcffilter --or -g "AO / RO > $FRACTION_AO_RO | AO / RO = $FRACTION_AO_RO" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" | sed -E -e 's#1\/0|0\/1|0\|1|1\|0#0|1#g' >> $TF/$SAMPLE.filtered.htz.heterozygous) \
>(vcffilter -g "RO / AO < $FRACTION_AO_RO" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' >> $TF/$SAMPLE.filtered.htz.alt-hmz) |\
vcffilter -g "AO / RO < $FRACTION_AO_RO" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#0\/1|1\/0|0\|1|1\|0|0\/2|2\/0|0\|2|2\|0#0|0#g' >> $TF/$SAMPLE.filtered.htz.ref-hmz

FALSEHTZALT=$(cat $TF/$SAMPLE.filtered.htz.alt-hmz | grep -v "^#" | wc -l) && echo "false HTZ (ALT): "$FALSEHTZALT
FALSEHTZREF=$(cat $TF/$SAMPLE.filtered.htz.ref-hmz | grep -v "^#" | wc -l) && echo "false HTZ (REF): "$FALSEHTZREF
TRUEHTZ=$(cat $TF/$SAMPLE.filtered.htz.heterozygous | grep -v "^#" | wc -l) && echo "true HTZ: "$TRUEHTZ
HTZSUM=$(echo "$FALSEHTZALT+$FALSEHTZREF+$TRUEHTZ" | bc -l)
echo "$HTZSUM of $HTZ processed/sorted"
#false HTZ (ALT): 29327
#false HTZ (REF): 37773
#true HTZ: 559952
#627052 of 627052 processed/sorted


#############################################################  not in the biallelic filtered dataset
#### NUMALT=2,3etc --> multiple ALT alleles

#### HTZ2

cat $TF/$SAMPLE.tmp.htz2 | tee \
>(vcffilter --or -g "RO / DP > $FRACTION_total | RO / DP = $FRACTION_total | RO / DP < $FRACTIONuLIMIT | RO / DP = $FRACTIONuLIMIT" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" >> $TF/$SAMPLE.filtered.htz2.heterozygous) \
>(vcffilter -g "RO / DP < $FRACTION_total" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' | sed -E -e 's#\/3|3\/0|0\|3|3\|0#3|3#g' >> $TF/$SAMPLE.filtered.htz2.alt-hmz) |\
vcffilter -g "RO / DP > $FRACTIONuLIMIT" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#0\/1|1\/0|0\|1|1\|0|0\/2|2\/0|0\|2|2\|0|0\/3|3\/0|0\|3|3\|0#0|0#g' >> $TF/$SAMPLE.filtered.htz2.ref-hmz

HTZ2REFHMZ=$(cat $TF/$SAMPLE.filtered.htz2.ref-hmz | grep -v "#" | wc -l) && echo "HTZ2 = REF HMZ: "$HTZ2REFHMZ
HTZ2ALTHMZ=$(cat $TF/$SAMPLE.filtered.htz2.alt-hmz | grep -v "#" | wc -l) && echo "HTZ2 = ALT HMZ: "$HTZ2ALTHMZ
HTZ2HTZ=$(cat $TF/$SAMPLE.filtered.htz2.heterozygous | grep -v "#" | wc -l) && echo "HTZ2 = HTZ: "$HTZ2HTZ
HTZ2SUM=$(echo "$HTZ2REFHMZ+$HTZ2ALTHMZ+$HTZ2HTZ" | bc -l)
echo "processed HTZ2 variants $HTZ2SUM of $HTZ2"


#### HTZ3

cat $TF/$SAMPLE.tmp.htz3 | tee \
>(vcffilter --or -g "RO / DP > $FRACTION_total | RO / DP = $FRACTION_total | RO / DP < $FRACTIONuLIMIT | RO / DP = $FRACTIONuLIMIT" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" >> $TF/$SAMPLE.filtered.htz3.heterozygous) \
>(vcffilter -g "RO / DP < $FRACTION_total" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#1\|2|1\|3|1\/2|1\/3#1|1#g' | sed -E -e 's#2\|1|2\|3|2\/1|2\/3#2|2#g' | sed -E -e 's#3\|1|3\|2|3\/1|3\/2#3|3#g' >> $TF/$SAMPLE.filtered.htz3.alt-hmz) |\
vcffilter -g "RO / DP > $FRACTIONuLIMIT" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#1\/2|2\/1|1\|2|2\|1|1\/3|3\/1|1\|3|3\|1|2\/3|3\/2|2\|3|3\|2#0|0#g' >> $TF/$SAMPLE.filtered.htz3.ref-hmz

HTZ3REFHMZ=$(cat $TF/$SAMPLE.filtered.htz3.ref-hmz | grep -v "#" | wc -l) && echo "HTZ3 = REF HMZ: "$HTZ3REFHMZ
HTZ3ALTHMZ=$(cat $TF/$SAMPLE.filtered.htz3.alt-hmz | grep -v "#" | wc -l) && echo "HTZ3 = ALT HMZ: "$HTZ3ALTHMZ
HTZ3HTZ=$(cat $TF/$SAMPLE.filtered.htz3.heterozygous | grep -v "#" | wc -l) && echo "HTZ3 = HTZ: "$HTZ3HTZ
HTZ3SUM=$(echo "$HTZ3REFHMZ+$HTZ3ALTHMZ+$HTZ3HTZ" | bc -l)
echo "processed HTZ3 variants $HTZ3SUM of $HTZ3"




#### HTZ partially mising    !! implements check/filter for too low cov of either REF or ALT

#cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -f10 | cut -d":" -f1 | sort | uniq
#.|0
#.|1

#cat $TF/$SAMPLE.tmp.partiallymissing |\
#grep -v "#" | cut -f10 | sort | uniq | tail
#
#.|1:6.93663e-05:5:0,1,0:0:0:1,0:37,0:-3.33503,-0.30103,0,-3.33503,-0.30103,-3.33503
#.|1:7.00617e-05:5:0,0,2:0:0:0,2:0,74:-6.65609,-6.65609,-6.65609,-0.60206,-0.60206,0
#.|1:7.04154e-12:16:3,9,0:3:107:9,0:315,0:-24.6492,0,-5.97156,-25.5523,-8.68083,-34.2392
#.|1:7.05518e-05:6:0,4:0:0:4:141:-12.6356,-1.80618,0
#.|1:7.36515e-06:6:1,3:1:37:3:111:-8.7554,0,-2.09398
#.|1:7.42988e-06:5:0,3,0,0:0:0:3,0,0:107,0,0:-9.61734,-0.90309,0,-9.61734,-0.90309,-9.61734,-9.61734,-0.90309,-9.61734,-9.61734
#.|1:7.46414e-05:7:0,5:0:0:5:185:-16.7115,-1.50515,0
#.|1:7.60118e-07:6:0,5:0:0:5:185:-16.6227,-1.50515,0

# fix multiple ALT genotypes, see above
# for now implemented based on RO / DP
## !! problematic is still how to handle ALT1/ALT2 cases, they should be rare, thus not implemented here

##  split vcf
#cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f1-14 > $TF/$SAMPLE.tmp.partiallymissing-a
#cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f15 > $TF/$SAMPLE.tmp.partiallymissing-b
#cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f16- > $TF/$SAMPLE.tmp.partiallymissing-c
##  for the AO field, take the highest AO count, parallel calc the sum of AO counts (all this is relevant for several ALT alleles only, frequently left over after decomposition)
### runs v slow, no alternatives yet found
#cat $TF/$SAMPLE.tmp.partiallymissing-b | parallel -k -j10 "echo {} | sed 's#,#\n#g' | tee >(sort -rn | head -n1 >> $TF/$SAMPLE.tmp.partiallymissing-b-highest) | numsum -x1 >> $TF/$SAMPLE.tmp.partiallymissing-b-sum"
##  paste together again
##paste $TF/$SAMPLE.tmp.partiallymissing-b $TF/$SAMPLE.tmp.partiallymissing-b-highest | head -n10000
#cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.partiallymissing-modified.highest
#cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.partiallymissing-modified.sum
#paste -d":" $TF/$SAMPLE.tmp.partiallymissing-a $TF/$SAMPLE.tmp.partiallymissing-b-highest $TF/$SAMPLE.tmp.partiallymissing-c >> $TF/$SAMPLE.tmp.partiallymissing-modified.highest
#paste -d":" $TF/$SAMPLE.tmp.partiallymissing-a $TF/$SAMPLE.tmp.partiallymissing-b-sum $TF/$SAMPLE.tmp.partiallymissing-c >> $TF/$SAMPLE.tmp.partiallymissing-modified.sum
#rm -f $TF/$SAMPLE.tmp.partiallymissing-a $TF/$SAMPLE.tmp.partiallymissing-b $TF/$SAMPLE.tmp.partiallymissing-c $TF/$SAMPLE.tmp.partiallymissing-b-highest $TF/$SAMPLE.tmp.partiallymissing-b-sum

#cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -f10 | cut -d":" -f1 | sort | uniq

cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f1-14 > $TF/$SAMPLE.tmp.partiallymissing-a
cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f15 > $TF/$SAMPLE.tmp.partiallymissing-b
cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "#" | cut -d":" -f16- > $TF/$SAMPLE.tmp.partiallymissing-c


cat $TF/$SAMPLE.tmp.partiallymissing-b | sed "s#,0##g" > $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes
cat $TF/$SAMPLE.tmp.partiallymissing-b  | sed "s#,0##g"  | awk -F ',' ' { print NF-1 } ' > $TF/$SAMPLE.tmp.partiallymissing.commanumbers
MAXCOMMA=$(cat $TF/$SAMPLE.tmp.partiallymissing.commanumbers | numbound) && echo $MAXCOMMA
cat $TF/$SAMPLE.tmp.partiallymissing.commanumbers | awk -v MAXCOMMA="$MAXCOMMA" ' { print MAXCOMMA-$1}' > $TF/$SAMPLE.tmp.partiallymissing.commas2add
cat $TF/$SAMPLE.tmp.partiallymissing.commas2add | sed "s#1#,#g" | sed "s#2#,,#g" | sed "s#3#,,,#g" | sed "s#4#,,,,#g" | sed "s#5#,,,,,#g" | sed "s#6#,,,,,,#g" | sed "s#7#,,,,,,,#g" | sed "s#8#,,,,,,,,#g" > $TF/$SAMPLE.tmp.partiallymissing.commas2substituted
paste -d "" $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes $TF/$SAMPLE.tmp.partiallymissing.commas2substituted > $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes2

Rscript $HOME/scripts/maxnumber.R $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes2
cat $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes2.csv | tail -n+2 | cut -d"," -f1 > $TF/$SAMPLE.tmp.partiallymissing.fixedcolumn

cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.partiallymissing-modified.highest
paste -d":" $TF/$SAMPLE.tmp.partiallymissing-a $TF/$SAMPLE.tmp.partiallymissing.fixedcolumn $TF/$SAMPLE.tmp.partiallymissing-c >> $TF/$SAMPLE.tmp.partiallymissing-modified.highest

rm -f $TF/$SAMPLE.tmp.partiallymissing.commanumbers $TF/$SAMPLE.tmp.partiallymissing.commas2add $TF/$SAMPLE.tmp.partiallymissing.commas2substituted $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes2 $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes
rm -f $TF/$SAMPLE.tmp.partiallymissing-a $TF/$SAMPLE.tmp.partiallymissing-b $TF/$SAMPLE.tmp.partiallymissing-c $TF/$SAMPLE.tmp.partiallymissing-b-highest $TF/$SAMPLE.tmp.partiallymissing.fixgenotypes2.csv




#cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | grep -v "#" | cut -f10 | cut -d":" -f1 | sort | uniq
rm -f $TF/$SAMPLE.filtered.partiallymissing.heterozygous $TF/$SAMPLE.filtered.partiallymissing.alt-hmz $TF/$SAMPLE.filtered.partiallymissing.ref-hmz $TF/$SAMPLE.filtered.partiallymissing.missing
cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | tee \
>(vcffilter --or -g "RO / AO > $FRACTION_AO_RO | RO / AO = $FRACTION_AO_RO" |\
vcffilter --or -g "AO / RO > $FRACTION_AO_RO | AO / RO = $FRACTION_AO_RO" |\
vcffilter -g "AO + RO > $COVEXCL" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -r 's#:GL\t\.\|0:|:GL\t\.\|1:#:GL\t0|1:#g' |\
sed -r 's#:GL\t\.\|2:#:GL\t0|2:#g' |\
sed -r 's#:GL\t\.\|3:#:GL\t0|3:#g' >> $TF/$SAMPLE.filtered.partiallymissing.heterozygous) \
>(vcffilter -g "RO / AO < $FRACTION_AO_RO" |\
vcffilter -g "AO + RO > $COVEXCL" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#\.\|0|0\|\.|\.\|1|1\|\.|\.\|2|2\|\.|\.\|3|3\|\.#1|1#g' >> $TF/$SAMPLE.filtered.partiallymissing.alt-hmz) \
>(vcffilter -g "AO + RO < $COVINCL" | grep -v "#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" | sed -r 's#:GL\t\.\|0:|:GL\t\.\|1:#:GL\t.:#g' >> $TF/$SAMPLE.filtered.partiallymissing.missing) |\
vcffilter -g "AO / RO < $FRACTION_AO_RO" |\
vcffilter -g "AO + RO > $COVEXCL" |\
grep -v "^#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" |\
sed -E -e 's#\.\|0|0\|\.|\.\|1|1\|\.|\.\|2|2\|\.|\.\|3|3\|\.#0|0#g' >> $TF/$SAMPLE.filtered.partiallymissing.ref-hmz


MISSINGGENOTYPES_partiallyMISSING=$(cat $TF/$SAMPLE.filtered.partiallymissing.missing | grep -v "#" | wc -l)
echo "additionally missing: "$MISSINGGENOTYPES_partiallyMISSING


PARTMISSREFHMZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.ref-hmz | grep -v "#" | wc -l) && echo "partially missing = REF HMZ: "$PARTMISSREFHMZ
PARTMISSALTHMZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.alt-hmz | grep -v "#" | wc -l) && echo "partially missing = ALT HMZ: "$PARTMISSALTHMZ
PARTMISSHTZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.heterozygous | grep -v "#" | wc -l) && echo "partially missing = HTZ: "$PARTMISSHTZ
PARTMISSINGSUM=$(echo "$PARTMISSREFHMZ+$PARTMISSALTHMZ+$PARTMISSHTZ+$MISSINGGENOTYPES_partiallyMISSING" | bc -l)
echo "processed partially missing variants $PARTMISSINGSUM of $PARTMISSING"


#partially missing = REF HMZ: 713
#partially missing = ALT HMZ: 680
#partially missing = HTZ: 102
#processed partially missing variants 1495 of 1495

#check if everything is caught
#cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | grep -v "#" | cut -f1,2 | sort -n > test1
#cat $TF/$SAMPLE.filtered.partiallymissing.ref-hmz | grep -v "#" | cut -f1,2 > test2
#cat $TF/$SAMPLE.filtered.partiallymissing.alt-hmz | grep -v "#" | cut -f1,2 >> test2
#cat $TF/$SAMPLE.filtered.partiallymissing.heterozygous | grep -v "#" | cut -f1,2 >> test2
#cat test2 | sort -n > test3

#diff test1 test3                                                                                                                       
#cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | vcffilter -g "AO + RO < $COVINCL" | grep -v "#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" | sed -r 's#:GL\t\.\|0:|:GL\t\.\|1:#:GL\t.:#g' | grep -A6 -P Group9.12"\t"669118
#cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | vcffilter -g "AO + RO > $COVEXCL" | grep -v "#" | grep -vP "GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$" 
#cat $TF/$SAMPLE.tmp.partiallymissing-modified.highest | grep -v "#" | grep -A6 -P Group9.12"\t"669118
#.|0:160.002:6:6,0:0:0:0:0:0,-2.66695,-4.47117
#.|0:160.002:6:6,0:0:0:0:0:0,-2.66695,-4.47117


#####
echo "heterozygous VCFs processed"



##################################################
##################################################
##################################################

echo "\
IN nonmissing missing sum \n\
$LINESBEFORE input\n\
$NONMISSINGGENOTYPES nonmissing\n\
$MISSINGGENOTYPES missing\n\
$FIRSTSUM sum\n\
 \n\
NONMISSINGGENOTYPES HMZ HTZ HTZ2/3(multialelic) PartiallyMissing SUM \n\
$HMZ homozygous\n\
$HTZ heterozygous\n\
$HTZ2 multiallelic REF/ALT\n\
$HTZ3 multiallelic ALT/ALT\n\
$PARTMISSING partially mising GT\n\
$SECONDSUM sum\n\
 \n\
HTZ REF-HMZ ALT_HMZ HTZ SUM \n\
$FALSEHTZREF heterozygous locus actually homozygous REF\n\
$FALSEHTZALT heterozygous locus actually homozygous ALT\n\
$TRUEHTZ true heterozygous\n\
$HTZSUM sum\n\
 \n\
HTZ2 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ2REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ2ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ2HTZ true heterozygous\n\
$HTZ2SUM sum\n\
 \n\
HTZ3 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ3REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ3ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ3HTZ true heterozygous\n\
$HTZ3SUM sum\n\
 \n\
partially missing REF-HMZ ALT_HMZ HTZ SUM \n\
$PARTMISSREFHMZ partially missing GT locus actually homozygous REF\n\
$PARTMISSALTHMZ partially missing GT locus actually homozygous ALT\n\
$PARTMISSHTZ partially missing GT locus actually heterozygous\n\
$MISSINGGENOTYPES_partiallyMISSING partially missing GT set to missing for too low coverage\n\
$PARTMISSINGSUM sum\n\
"

echo "\
IN nonmissing missing sum \n\
$LINESBEFORE input\n\
$NONMISSINGGENOTYPES nonmissing\n\
$MISSINGGENOTYPES missing\n\
$FIRSTSUM sum\n\
 \n\
NONMISSINGGENOTYPES HMZ HTZ HTZ2/3(multialelic) PartiallyMissing SUM \n\
$HMZ homozygous\n\
$HTZ heterozygous\n\
$HTZ2 multiallelic REF/ALT\n\
$HTZ3 multiallelic ALT/ALT\n\
$PARTMISSING partially mising GT\n\
$SECONDSUM sum\n\
 \n\
HTZ REF-HMZ ALT_HMZ HTZ SUM \n\
$FALSEHTZREF heterozygous locus actually homozygous REF\n\
$FALSEHTZALT heterozygous locus actually homozygous ALT\n\
$TRUEHTZ true heterozygous\n\
$HTZSUM sum\n\
 \n\
HTZ2 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ2REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ2ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ2HTZ true heterozygous\n\
$HTZ2SUM sum\n\
 \n\
HTZ3 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ3REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ3ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ3HTZ true heterozygous\n\
$HTZ3SUM sum\n\
 \n\
partially missing REF-HMZ ALT_HMZ HTZ SUM \n\
$PARTMISSREFHMZ partially missing GT locus actually homozygous REF\n\
$PARTMISSALTHMZ partially missing GT locus actually homozygous ALT\n\
$PARTMISSHTZ partially missing GT locus actually heterozygous\n\
$MISSINGGENOTYPES_partiallyMISSING partially missing GT set to missing for too low coverage\n\
$PARTMISSINGSUM sum\n\
" > $FILTEROUTFOLDER/$SAMPLE.filtered.txt

TOTALMISSING=$(echo "$MISSINGGENOTYPES+$MISSINGGENOTYPES_partiallyMISSING" | bc -l)
TOTALHTZ=$(echo "$PARTMISSHTZ+$HTZ2HTZ+$HTZ3HTZ+$TRUEHTZ" | bc -l)
TOTALHMZ=$(echo "$HMZ+$HTZ2REFHMZ+$HTZ2ALTHMZ+$HTZ3REFHMZ+$HTZ3ALTHMZ+$FALSEHTZREF+$FALSEHTZALT+$PARTMISSREFHMZ+$PARTMISSALTHMZ" | bc -l)
TOTAL=$(echo "$TOTALHTZ+$TOTALHMZ+$TOTALMISSING" | bc -l)

echo "total missing $TOTALMISSING"
echo "total HTZ $TOTALHTZ"
echo "total HMZ $TOTALHMZ"
echo "total $TOTAL from input $LINESBEFORE"


echo "total missing $TOTALMISSING" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total HTZ $TOTALHTZ" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total HMZ $TOTALHMZ" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total $TOTAL from input $LINESBEFORE" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt


##################################################

## combine pieces and rudimentary sort
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.filtered.combined
cat \
$TF/$SAMPLE.filtered.homozygous \
$TF/$SAMPLE.filtered.htz2.alt-hmz \
$TF/$SAMPLE.filtered.htz2.heterozygous \
$TF/$SAMPLE.filtered.htz2.ref-hmz \
$TF/$SAMPLE.filtered.htz3.alt-hmz \
$TF/$SAMPLE.filtered.htz3.heterozygous \
$TF/$SAMPLE.filtered.htz3.ref-hmz \
$TF/$SAMPLE.filtered.htz.alt-hmz \
$TF/$SAMPLE.filtered.htz.heterozygous \
$TF/$SAMPLE.filtered.htz.ref-hmz \
$TF/$SAMPLE.filtered.missing \
$TF/$SAMPLE.filtered.partiallymissing.alt-hmz \
$TF/$SAMPLE.filtered.partiallymissing.heterozygous \
$TF/$SAMPLE.filtered.partiallymissing.missing \
$TF/$SAMPLE.filtered.partiallymissing.ref-hmz |\
grep -v "#" | sed -r "s#\t\.\|\.#\t.#g" | sort -n >> $TF/$SAMPLE.filtered.combined

echo "combined VCF parts"
#cat $TF/$SAMPLE.filtered.combined | grep -v "#" | cut -f10 | cut -d":" -f1 | sort | uniq


# sort, format missing genotypes, set heterozygous to unphased
vt sort $TF/$SAMPLE.filtered.combined |\
vcffixup - |\
sed -r "s#\t0\|1:#\t0/1:#g" |\
sed -r "s#\t0\|2:#\t0/2:#g" | sed -r "s#\t0\|3:#\t0/3:#g" |\
sed -r "s#\t1\|2:#\t1/2:#g"| sed -r "s#\t2\|3:#\t2/3:#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz

FUSED=$(cat $TF/$SAMPLE.filtered.combined | grep -v "#" | wc -l)
FINAL=$(zcat $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz | grep -v "#" | wc -l)
echo "original: $LINESBEFORE ; fused: $FUSED ; final: $FINAL"
echo "original: $LINESBEFORE ; fused: $FUSED ; final: $FINAL" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt

#vt peek $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz
#stats: no. of samples                     :          1
#       no. of chromosomes                 :        214
#       no. of SNP                         :    2734578

#vt peek $VCFINPUT
#stats: no. of samples                     :          1
#       no. of chromosomes                 :        214
#       no. of SNP                         :    2734578

#cat $TF/$SAMPLE.filtered.header > temp
#zcat $VCFINPUT | grep -v "#" | sort -n >> temp
#vt sort temp > temp2
##  zcat $TF/$SAMPLE.filtered.vcf.gz | grep -v "#" | cut -f 1,2 > output
##  zcat $VCFINPUT | grep -v "#" | cut -f 1,2 > input
#cat temp2 | grep -v "#" | cut -f 1,2 > input2
#diff input output
##differences in sorting
#diff input2 output
## after resorting the input the same way, there are no differences between in/output in terms of scf/position/order


MISSINGBEFORE=$(zcat $VCFINPUT | grep -v "#" | grep -c "\.:\.:\.:\.:\.") && echo $MISSINGBEFORE
## MISSINGGTAFTER=$(zcat $FILF/$SAMPLE.filtered.vcf.gz | grep -v "#" | grep -c -e "\.:*:*:") && echo $MISSINGGTAFTER
MISSINGAFTER=$TOTALMISSING

echo "number of missing genotypes before/after/difference processing/filtering/setting Htz to missing"
MISSINGDIFF=$(echo "$MISSINGGENOTYPES-$MISSINGBEFORE" | bc -l)

#echo "SAMPLE : TOTAL : MISSINGGTBEFORE : MISSINGGTAFTER : MISSINGDIFF : TOTALHTZ : TOTALHMZ : HTZ : TRUEHTZ : FALSEHTZREF : FALSEHTZALT : HTZ2 : HTZ3"
#echo "$SAMPLE : $TOTAL : $MISSINGGTBEFORE : $MISSINGGTAFTER : $MISSINGDIFF : $TOTALHTZ : $TOTALHMZ : $HTZ : $TRUEHTZ : $FALSEHTZREF : $FALSEHTZALT : $HTZ2 : $HTZ3"
echo "$SAMPLE : $TOTAL : $MISSINGBEFORE : $MISSINGAFTER : $MISSINGDIFF : $TOTALHTZ : $TOTALHMZ : $HTZ : $TRUEHTZ : $FALSEHTZREF : $FALSEHTZALT : $HTZ2 : $HTZ3"\
>> $FILTEROUTFOLDER/stats.filtered.txt

#cleanup
cd $ORIGFOLDER
rm -rf $TF/

   #echo "start: $STARTTIME"
   #ENDTIME=$(date) && echo "end $ENDTIME"

echo "############################################################################"
echo "##########       $SAMPLE done"
echo "############################################################################  \n \n \n"
### END OF GT-SCRIPT
