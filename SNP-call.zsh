#!/bin/zsh
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

#CPUs=30
#INPUTFOLDER="$HOME/scratch/scratchspace/melipona/bam"
#OUTPUTFOLDER="$HOME/scratch/scratchspace/melipona/vcf_mqua"
#REF="$HOME/ref/ref_w_controls/Mqua.fa"
#cat $HOME/ref/Mqua/Mqua.fa $HOME/ref/Mqua/NC_026198_Mscut_mt.fasta $HOME/ref/wolbachia/wolbachia.fasta $HOME/ref/phiX/NC_001422.fasta > $HOME/ref/ref_w_controls/Mqua.fa
#bwa index $HOME/ref/ref_w_controls/Mqua.fa
#$HOME/scripts/mapping_bwa.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs $REF

#exclude mitochondrion
#MT NC_026198.1 Mscut
#MT="NC_026198"
#NC_001566.1 Apis_mellifera_ligustica_mitochondrion
#CP001391.1 Wolbachia_sp_wRi
#AM999887.1 Wolbachia_Culex_quinquefasciatus_Pel_strain_wPip
#NC_001422.1 phiX

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


vt validate  $OUTPUTFOLDER/scaffolds_htz_fused/Mqua.GT.vcf.gz
#1 773 845



cat Q01.depth.HiCov.bed | grep KQ435687.1'\t'299473











#| vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED | vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | bgzip -@ 2 -f -c /dev/stdin > $OUTPUTFOLDER/scaffolds/{}.joint.vcf.gz && \
#tabix -fp vcf $OUTPUTFOLDER/scaffolds/{}.joint.vcf.gz"

#vt cat -n $OUTPUTFOLDER/scaffolds/*.joint.vcf.gz -o $OUTPUTFOLDER/combine.joint.vcf
#bgzip -f -@ $CPUs -c $OUTPUTFOLDER/combine.joint.vcf > $OUTPUTFOLDER/combine.joint.vcf.gz && tabix -fp vcf $OUTPUTFOLDER/combine.joint.vcf.gz
#vt peek $OUTPUTFOLDER/combine.joint.vcf.gz > $OUTPUTFOLDER/combine.stats

#zcat ${OUTPUTFOLDER}/scaffolds_SSRfiltered/KQ435840.1:0-3151304.vcf.gz | vcffilter -f '! ( TYPE = ins | TYPE = del | TYPE = complex )' |\
#vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
#grep -v -P 'TYPE=ins|TYPE=del|TYPE=complex' |\
#vcfstreamsort -a | vcfuniq | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_filtered/KQ435840.1:0-3151304.snps.vcf.gz && \
#tabix -fp vcf $OUTPUTFOLDER/scaffolds_filtered/KQ435840.1:0-3151304.snps.vcf.gz

#zcat $OUTPUTFOLDER/scaffolds_filtered/*.snps.vcf.gz | vcffirstheader | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/Mqua.tmp.vcf.gz

zcat $OUTPUTFOLDER/Mqua.tmp.vcf.gz | vcfstreamsort | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/Mqua.snps.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/Mqua.snps.vcf.gz

rm -f $OUTPUTFOLDER/Mqua.snps.vcf.gz

vcfstats $OUTPUTFOLDER/Mqua.snps.vcf.gz
vt peek $OUTPUTFOLDER/Mqua.snps.vcf.gz
#stats: no. of samples                     :         53
#       no. of chromosomes                 :        713
#       no. of SNP                         :   10216499
#           2 alleles                      :        10136906 (2.94) [7567142/2569764]
#           3 alleles                      :           78748 (0.55) [55685/101811]
#           4 alleles                      :             845 (0.50) [844/1691]


dumpContigsFromHeader
popStats --type GL --target 0,1,2,3,4,5,6,7 --file
#FATAL: unknown genotype: .|.
vcfcountalleles
vcfhetcount
vcfaltcount
vcfhethomratio
vcfinfosummarize
vcfnumalt ##stream # NUMALT based on sample GT
vcfnull2ref
vcfglxgt
vcfnulldotslashdot
vcfgenoexpandnull https://github.com/vcflib/vcflib/pull/92
vcfstats
vcfsitesummarize KQ436020.1:0-6715.snps.vcf.gz

zcat KQ438543.1:0-1874.vcf.gz | vcfnull2ref - | vcfnumalt -
zcat KQ436014.1:0-4520.snps.vcf.gz | grep -v "#" | head -n100
popStats --type GL --target 0,1,2,3,4,5,6,7 --file KQ436014.1:0-4520.snps.vcf.gz

zcat KQ436014.1:0-4520.snps.vcf.gz | vcffilter -f '( AB > 0 )'


mkdir -p $OUTPUTFOLDER/scaffolds_speciesdiff
mkdir -p $OUTPUTFOLDER/scaffolds_htz
cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_filtered/{}.snps.vcf.gz | vcffilter -f '( AB = 0 )' |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_speciesdiff/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_speciesdiff/{}.vcf.gz"
/home/es/scratch/scratchspace/melipona/vcf_mqua/scaffolds_speciesdiff/KQ435744.1:0-4956068.vcf.gz
Unsorted positions on sequence #1: 3578257 followed by 3578254
KQ435793.1:0-1549515
parallel: Warning: No more file handles. Raising ulimit -n or /etc/security/limits.conf may help.

cat $TARGETS | parallel --no-notice -k -j $CPUs \
"echo {} && zcat ${OUTPUTFOLDER}/scaffolds_filtered/{}.snps.vcf.gz | vcffilter -f '( AB > 0 )' |\
bgzip -f -@ 3 -c /dev/stdin > $OUTPUTFOLDER/scaffolds_htz/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/scaffolds_htz/{}.vcf.gz"
[E::hts_idx_push] Unsorted positions on sequence #1: 1259531 followed by 1259523
tbx_index_build failed: /home/es/scratch/scratchspace/melipona/vcf_mqua/scaffolds_htz/KQ435714.1:0-1906513.vcf.gz




vcfstats $OUTPUTFOLDER/scaffolds_speciesdiff/KQ438551.1:0-461204.vcf.gz
#11196
vcfstats $OUTPUTFOLDER/scaffolds_htz/KQ438551.1:0-461204.vcf.gz
#total variant sites:	2574 of which 2518 (0.978244) are biallelic and 56 (0.021756) are multiallelic, total variant alleles:	2630

vcfglxgt KQ438543.1:0-1874.vcf.gz
zcat KQ438543.1:0-1874.vcf.gz | vcfglxgt
zcat KQ438543.1:0-1874.vcf.gz | vcfnull2ref - | vcfnumalt -
zcat KQ438543.1:0-1874.vcf.gz | vcfnumalt -
0:.:.:.:.:.:.:.
0/1:56:48,8:48:1656:8:261:-6.95035,0,-132.409
28%
0/1:	38:	8,20,7:	8:	256:	20:	500:	-36.7182,-0.233497,-14.809
GT:	DP:	AD:	RO:	QR:	AO:	QA:	GL
zcat KQ438543.1:0-1874.vcf.gz | vcfnull2ref - | vcfnumalt - | vcffilter -g '( ( DP ) > 4 )' | vcffilter -t -a -k -g '! ( ( RO / DP ) < 0.35 )'
#vcffilter -t -a -k -g '! ( ( AO / DP ) < 0.35 )'
| vcffilter -t -a -k -g '! ( ( RO / DP ) < 0.35 )'


### separate monomorphic sites (species differences)


## concat

##extract per sample
#rm -rf $PROJECT/mt/vcf_individuals_biallelicQ20
mkdir -p $PROJECT/mt/vcf_individuals_biallelicQ20
vcfsamplenames $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz {} | vcffixup - | grep -v -P AO:QA:GL'\t'0: | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals_biallelicQ20/MT.biallelicQ20.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals_biallelicQ20/MT.biallelicQ20.{}.vcf.gz"






## re-decompose
#zcat mastervcf/combine.joint.pops_no_dups.decomposed.noSSR.noGT.noInDel.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED | vcffixup - | #vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles > mastervcf/tmp.vcf
#bgzip -@ $CPUs -f -c mastervcf/tmp.vcf > mastervcf/combine.joint.pops_no_dups.decomposed.noSSR.noGT.noInDel.decomposed.vcf.gz && \
#tabix -fp vcf mastervcf/combine.joint.pops_no_dups.decomposed.noSSR.noGT.noInDel.decomposed.vcf.gz
#rm -f mastervcf/tmp.vcf

zcat $OUTPUTFOLDER/combine.filtered2.vcf.gz | grep -v "#" | cut -f4,5 | sort | uniq

vt validate $OUTPUTFOLDER/combine.filtered2.vcf.gz
vt validate $OUTPUTFOLDER/combine.InDel.vcf.gz

rm -f $OUTPUTFOLDER/combine.joint.vcf

vt peek $OUTPUTFOLDER/combine.filtered2.vcf.gz

#####################




















