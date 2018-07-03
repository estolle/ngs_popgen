#!/bin/zsh

# incomplete Script at the end


############################################################################
# SNPcall with freebayes + R viz
# Eckart Stolle, Sept 2017
# limit for SNP call set to 1000 bp scaffolds (change if not appopriate)
############################################################################
# uses zsh shell (not bash), get it by: sudo apt-get install zsh
# usage: $HOME/scripts/SNP-call.zsh INPUTFOLDER OUTPUTFOLDER CPUs REFERENCE


if [ $# -ne 4 ]; then
    echo $0: usage: ./SNP-call.zsh INPUTFOLDER OUTPUTFOLDER CPUs REFERENCE
	echo "\nINPUTFOLDER"
	echo "\nOUTPUTFOLDER"
	echo "\nCPUs: threads/cores to be used"
	echo "\nREFERENCE: reference (genome fasta) to be used"
    exit 1
fi

#set/get variables
INPUTFOLDER=$1
OUTPUTFOLDER=$2
CPUs=$3
REF=$4

MT="NC_001566.1"
PHIX="NC_001422.1"


# limit scafffold length
LENGTHLIMIT=1000

## SNP call thresholds
MINALTFRAC=0.35
#MINALTFRAC=0.20
MINALTN=4
MINCOV=8



####uses these progz:
which freebayes
which vcffilter
which vcfallelicprimitives
which vcffixup
which vcftools
which vt
which tabix
which bgzip
which R


#if some inputfiles do not exist
if [[ ! -f $(ls -1 $INPUTFOLDER/*.bam | head -n1) ]]; then
echo "input file(s) does not exist"
    exit 1
fi
if [[ -f $(ls -1 $INPUTFOLDER/*.bam | head -n1) ]]; then
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


echo "SNP calling"

# make outputfolder for indiv. scaffolds
mkdir -p $OUTPUTFOLDER/scaffolds_unfiltered/

# fetch list of files to process
ls -1 $INPUTFOLDER/*.bam | rev | cut -d"." -f 2- | cut -d"/" -f 1 | rev > $OUTPUTFOLDER/samples.lst
ls -1 $INPUTFOLDER/*.bam > $OUTPUTFOLDER/variantsamples.lst

LISTtoPROCESS="$OUTPUTFOLDER/variantsamples.lst"

#limit scafffolds to size threshold, exclude Mitochondrion
cat $REF.region | grep -v $MT | tr '-' '\t' | awk -v LENGTHLIMIT="$LENGTHLIMIT" '$2 > LENGTHLIMIT' | tr '\t' '-' > $REF.region$LENGTHLIMIT


N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
echo "calling $N samples (bam files)"



#dryrun
echo "settings: \n xxxxx \n xxxx \n xxxx \n xxxx \n xxx \n xxxx \n xxxxxx"
sleep 3s
cat $LISTtoPROCESS | parallel --no-notice -j 1 "echo $INPUTFOLDER/{}.bam"
sleep 3s

######################### limit SNP call on region larger than 1000 bp

TARGETS="$REF.region$LENGTHLIMIT"
BAMLST=$LISTtoPROCESS
#MINALTFRAC=0.20
#MINALTFRAC=0.35
#MINALTN=4
#MINCOV=8
# alternativ run method freebayes-parallel $TARGETS $CPUs \
# filtering / threshold can be relaxed, as the MasterSNP list is used for GT and filtered there

cat $TARGETS | parallel --no-notice -k -j $CPUs \
"freebayes \
--region {} \
--fasta-reference $REF \
--bam-list $BAMLST \
--ploidy 2 \
--use-best-n-alleles 4 \
--use-reference-allele \
--haplotype-length -1 \
--max-complex-gap -1 \
--min-mapping-quality 45 \
--min-base-quality 20 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV | bgzip -@ 2 -f -c /dev/stdin > $OUTPUTFOLDER/scaffolds_unfiltered/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_unfiltered/{}.vcf.gz"


CPUs=12
mkdir -p $OUTPUTFOLDER/scaffolds_SSRfiltered
SSRBED="$REF.Mqua.fa.SSRs.ext5bp.bed"
HICOVBED="depth.all.cov.25.HiCov.bed"
#87297 bp
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_unfiltered/{}.vcf.gz | vcfintersect -v -l -b $SSRBED | vcfintersect -v -l -b $HICOVBED | vcffilter -f 'QUAL > 30' |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_SSRfiltered/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_SSRfiltered/{}.vcf.gz"

mkdir -p $OUTPUTFOLDER/scaffolds_indels
#grep -v -P 'TYPE=snp|TYPE=mnp' | grep -v 'TYPE=complex'
#vcffilter -f '( TYPE = ins | TYPE = del | TYPE = complex )'
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_SSRfiltered/{}.vcf.gz | vcfindels |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_indels/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_indels/{}.vcf.gz"


mkdir -p $OUTPUTFOLDER/scaffolds_snps
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_SSRfiltered/{}.vcf.gz | vcfnoindels | vcffilter -f '! ( TYPE = complex )'|\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED | vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - |\
vcfuniqalleles | vcfnoindels | vcfbiallelic | vcfnumalt - | vcfnulldotslashdot | vcfuniq |\
grep -v -P 'TYPE=ins|TYPE=del|TYPE=complex' | bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_snps/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_snps/{}.vcf.gz"
#sed -r 's/(.)\/(.)/\1|\2/g'


mkdir -p $OUTPUTFOLDER/scaffolds_speciesdiff
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_snps/{}.vcf.gz | vcffilter -f '( AB = 0 )' |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_speciesdiff/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_speciesdiff/{}.vcf.gz"


mkdir -p $OUTPUTFOLDER/scaffolds_htz
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_snps/{}.vcf.gz | vcffilter -f '( AB > 0 )' |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_htz/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_htz/{}.vcf.gz"

mkdir -p $OUTPUTFOLDER/scaffolds_htz_fused
zcat $OUTPUTFOLDER/scaffolds_htz/*.vcf.gz | vcffirstheader | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_htz_fused/Mqua.GT.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/scaffolds_htz_fused/Mqua.GT.vcf.gz












