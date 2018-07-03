#!/bin/zsh
############################################################################
# genotyping bam's based on SNPs in a master-SNP
# Eckart Stolle, Sept 2017, May 2018
############################################################################
################################################################################
#####  use a VCF to Genotype bam's

# uses zsh shell (not bash), get it by: sudo apt-get install zsh
# usage: ./SNP-genotyping.zsh INPUTFOLDER OUTPUTFOLDER CPUs REF VCFMASTER


#Step1 genotyping + basic filtering/decomposition of haplotypes, only biallelic are retained in this version

MINALTFRAC="0.25"
MINALTN=2
MINCOV=2


if [ $# -ne 5 ]; then
    echo $0: usage: ./SNP-genotyping.zsh INPUTFOLDER OUTPUTFOLDER CPUs REF VCFMASTER 
	echo "\nINPUTFOLDER: Folder containing bam (+bai) files"
	echo "\nOUTPUTFOLDER: Folder for genotype-vcfs , will be created new if not existing"
	echo "\nCPUs: threads/cores to be used"
	echo "\nREF: reference fasta"
	echo "\nVCF: master-vcf"
    exit 1
fi


#set/get variables
INPUTFOLDER=$1
OUTPUTFOLDER=$2
CPUs=$3
REF=$4
VCF=$5


#INPUTFOLDER="/scratch/ek/apis/thelytoky/bam/clones"
#OUTPUTFOLDER="/scratch/ek/apis/thelytoky/vcf/clones"
#CPUs=100
#REF="/scratch/genomes/genomes_with_controls_MT_phiX_Wolbachia/Amel_4.5_scaffolds.fa"
###master vcf
#VCF="$HOME/scratch/apis/thelytoky/vcf/4_merged_vcf/selection.vcfcombine.biallelic.vcf.gz"

#####
FOLDERNAME=$(echo $INPUTFOLDER | rev | cut -d"/" -f 1 | rev) && echo $FOLDERNAME
FILTEROUTFOLDER="$OUTPUTFOLDER/genotyping/filtered.sample"
#FOLDERNAME="clones"

mkdir -p $OUTPUTFOLDER
mkdir -p $OUTPUTFOLDER/genotyping/raw
mkdir -p $OUTPUTFOLDER/genotyping/raw.filtered
mkdir -p $OUTPUTFOLDER/genotyping/raw.filtered.fused
mkdir -p $OUTPUTFOLDER/genotyping/raw.filtered.fused.sample

echo "CPUs" $CPUs
CPUs23=$(printf "%.0f" $(echo "scale=2;($CPUs/3)*2" | bc)) && echo $CPUs23
CPUs2=$(printf "%.0f" $(echo "scale=2;$CPUs/2" | bc)) && echo $CPUs2
CPUs3=$(printf "%.0f" $(echo "scale=2;$CPUs/3" | bc)) && echo $CPUs3
CPUs4=$(printf "%.0f" $(echo "scale=2;$CPUs/4" | bc)) && echo $CPUs4
CPUs5=$(printf "%.0f" $(echo "scale=2;$CPUs/5" | bc)) && echo $CPUs5
CPUs6=$(printf "%.0f" $(echo "scale=2;$CPUs/8" | bc)) && echo $CPUs6

## make file of regions (based on .fai of the Ref)
cat $REF.fai | cut -f 1,2 | sed "s/\t/:0-/g" > $OUTPUTFOLDER/regions

# make lists of files to process
ls -1 $INPUTFOLDER/*.bam | rev | cut -d"." -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/$FOLDERNAME.samples.lst
ls -1 $INPUTFOLDER/*.bam > $OUTPUTFOLDER/$FOLDERNAME.bams.lst


LISTtoPROCESS="$OUTPUTFOLDER/$FOLDERNAME.samples.lst"
N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
echo "processing $N samples"
BAMLST="$OUTPUTFOLDER/$FOLDERNAME.bams.lst"


## prep the input vcf (fix so that its not done every time)
cat $OUTPUTFOLDER/regions | sed "s/:0-/\t0\t/g" > $OUTPUTFOLDER/bed

if [[ -d "$VCF.master.scf" ]]; then
echo "per-scaffold-vcfs for Reference VCF exist already, proceeding to next step"
else 
	echo "per-scaffold-vcfs for Reference VCF does not exist, ... creating"
    mkdir -p $VCF.master.scf 
	cat $OUTPUTFOLDER/regions | grep -vP 'NC_001566|CP001391|AM999887|NC_001422|GroupUn' | parallel -k -j $CPUs "\
echo {}; echo {} > $VCF.master.scf/{}.region"
	cat $OUTPUTFOLDER/regions | grep -vP 'NC_001566|CP001391|AM999887|NC_001422|GroupUn' | parallel -k -j $CPUs "\
echo {}; echo {} | sed 's/:0-/\t0\t/g' > $VCF.master.scf/{}.bed"
	cat $OUTPUTFOLDER/regions | grep -vP 'NC_001566|CP001391|AM999887|NC_001422|GroupUn' | parallel -k -j $CPUs2 "\
echo {}; tabix -hR $VCF.master.scf/{}.bed $VCF | bgzip -f -@ 10 -c /dev/stdin > $VCF.master.scf/{}.vcf.gz; tabix -fp vcf $VCF.master.scf/{}.vcf.gz"
	ls $VCF.master.scf/*.vcf.gz
	echo "per-scaffold-vcfs for Reference VCF created"
fi


## genotype

### MODIFY min base quality if necessary
### MODIFY min mapping qual: Illumina HiSeq PE 150 bp, bwa mapping: 45-60; SOLiD bowtie mapping 250-255, SE 50bp bwa: 36-27, SE100 bp 45-60
# Cridland: 32
# new bam index for RNAseq mappings:
# ls -1 *.bam | parallel -k -j20 "echo {} && sambamba index -p -t 20 {} {}.bai; echo {}"
# Wallberg14 data: ERROR(freebayes): Couldn't find read group id (@RG tag) for BAM Alignment SRR1152244.3116971 at position 81 in sequence


ulimit -n
#1024
#each thread consumes (n bamfiles * 2)+3 file descriptors, 60*20=1300ish
ulimit -n 10240
ulimit -n

#ls -1 *.bam | parallel -k -j20 "echo {} && sambamba index -p -t 20 {} {}.bai; echo {}"
#Group4.15:0-1083 in VCF header but no SNP, causing empty file
cat $OUTPUTFOLDER/regions | grep -vP 'NC_001566|CP001391|AM999887|NC_001422|GroupUn|Group4\.15' | parallel --no-notice -j $CPUs2 "\
echo {}; freebayes \
--fasta-reference $REF \
--haplotype-basis-alleles $VCF.master.scf/{}.vcf.gz \
--region {} \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 8 \
--haplotype-length 0 \
--min-mapping-quality 45 \
--min-base-quality 20 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV \
--use-reference-allele \
--bam-list $BAMLST |\
bgzip -f -@ 2 -c /dev/stdin > $OUTPUTFOLDER/genotyping/raw/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/genotyping/raw/{}.vcf.gz && ls $OUTPUTFOLDER/genotyping/raw/{}.vcf.gz"

#empty vcf's, Group4.15 without header
# 4.15 also makes freebayes error
# VCF input: ##contig=<ID=Group4.15,length=1083> in header but not any genotypes
# zcat $VCF | grep "Group4\.15"
# zcat $VCF | grep "Group3\.17"

#-rw-rw-r-- 1 ek ek   72 Feb  1 15:09 Group12.12:0-962.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:10 Group13.13:0-1409.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:12 Group15.1:0-1130.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:16 Group3.17:0-1698.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:17 Group4.15:0-1083.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:22 Group6.19:0-1436.vcf.gz.tbi
#-rw-rw-r-- 1 ek ek   72 Feb  1 15:23 Group6.33:0-2665.vcf.gz.tbi

# decompose, filter, only keep biallelic
##!! sed -r 's/(.)\/(.)/\1|\2/g' sets everything to phased

cat $OUTPUTFOLDER/regions | grep -vP 'NC_001566|CP001391|AM999887|NC_001422|GroupUn|Group4\.15' | parallel --no-notice -j $CPUs3 \
"echo {}; ls $OUTPUTFOLDER/genotyping/raw/{}.vcf.gz; zcat $OUTPUTFOLDER/genotyping/raw/{}.vcf.gz | vcffilter -f 'QUAL > 30' |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | vcfnoindels | vcfbiallelic | vcfnumalt - |\
vcfnulldotslashdot | sed -r 's/(.)\/(.)/\1|\2/g' |\
bgzip -f -@ 4 -c /dev/stdin > $OUTPUTFOLDER/genotyping/raw.filtered/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/genotyping/raw.filtered/{}.vcf.gz; ls $OUTPUTFOLDER/genotyping/raw.filtered/{}.vcf.gz"


rm -f $OUTPUTFOLDER/genotyping/raw.filtered/Group4.15:0-1083.vcf.gz $OUTPUTFOLDER/genotyping/raw.filtered/Group4.15:0-1083.vcf.gz.tbi
#rm -f */genotyping/raw/Group4.15:0-1083.vcf.gz */genotyping/raw/Group4.15:0-1083.vcf.gz.tbi
#rm -f */genotyping/raw.filtered/Group4.15:0-1083.vcf.gz */genotyping/raw.filtered/Group4.15:0-1083.vcf.gz.tbi


# fused VCF per scaffold/chrom into one vcf
vcfcat $OUTPUTFOLDER/genotyping/raw.filtered/*.vcf.gz | vcffirstheader | bgzip -f -@ $CPUs3 -c /dev/stdin > $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz

#zcat $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz | grep -v "#" | cut -f1 | sort | uniq


# split by sample
vcfsamplenames $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz > $OUTPUTFOLDER/samplenames.raw.filtered.VCF.txt

cat $OUTPUTFOLDER/samplenames.raw.filtered.VCF.txt | parallel -j $CPUs4 \
"echo {}; vcfkeepsamples $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz {} | bgzip -f -@ 5 -c /dev/stdin > $OUTPUTFOLDER/genotyping/raw.filtered.fused.sample/$FOLDERNAME.{}.raw.filtered.vcf.gz &&\
tabix -fp vcf $OUTPUTFOLDER/genotyping/raw.filtered.fused.sample/$FOLDERNAME.{}.raw.filtered.vcf.gz; echo {}"

## fix output
vt peek $VCF &> $OUTPUTFOLDER/stats.$FOLDERNAME.input.VCF.txt
vt peek $OUTPUTFOLDER/genotyping/raw.filtered.fused/$FOLDERNAME.raw.filtered.fused.vcf.gz &> $OUTPUTFOLDER/stats.$FOLDERNAME.raw.filtered.VCF.txt
echo "done"

###### end of script, now filtering via script

