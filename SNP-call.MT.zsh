
CPUs=15

REHEADERSCRIPT="$HOME/scripts/vcfreheader.zsh"
REHEADERSCRIPT="$PROJECT/bin/vcfreheader.zsh"
#Mqua
#MT="NC_026198"

#Soli
MT="NC_014672"

cat $REF.fai | grep $MT > $REF.mt.fai
fasta_generate_regions.py $REF.mt.fai 500 > $REF.mt.500bp.regions

mkdir -p mt/vcf_master
ls -1 $PROJECT/bam/*bam > $PROJECT/mt/vcf_master/bam.lst

BAMLST="$PROJECT/mt/vcf_master/bam.lst"
MINALTFRAC=0.50
MINALTN=5
MINCOV=5

(cat $REF.mt.500bp.regions | parallel -k -j "$CPUs" freebayes --region {} \
--fasta-reference $REF \
--ploidy 1 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 20 \
--haplotype-length 0 \
--min-mapping-quality 40 \
--min-base-quality 20 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV \
--use-reference-allele \
--bam-list $BAMLST > $PROJECT/mt/vcf_master/mt_1.vcf
)
cat $PROJECT/mt/vcf_master/mt_1.vcf | vcffirstheader | vcfstreamsort -a |\
vcfuniq | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED | vcfstreamsort -a |\
bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/mt_2.vcf.gz && \
tabix -fp vcf $PROJECT/mt/vcf_master/mt_2.vcf.gz

bgzip -f -@ $CPUs $PROJECT/mt/vcf_master/mt_1.vcf
tabix -fp vcf $PROJECT/mt/vcf_master/mt_1.vcf.gz

zcat mt_2.vcf.gz | grep -v "#" | cut -f1-5 > $PROJECT/mt/vcf_master/mt_2.genotypes
#zgrep -n "#CHROM" -m 1 $PROJECT/mt/vcf_master/mt.vcf.gz | cut -d : -f 1
vcfsamplenames $PROJECT/mt/vcf_master/mt_2.vcf.gz > $PROJECT/mt/vcf_master/mt_2.samplenames

#reheader
$REHEADERSCRIPT $PROJECT/mt/vcf_master/mt_2.vcf.gz


SNPS=$(zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | cut -f1-5 | wc -l)
echo "SNPs $SNPS"
HMZ=$(zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | grep -v -P '0[/|]1' | cut -f1-5 | wc -l)
echo "SNPs without a heterozygous individual: $HMZ"
#3097 vs  4350
#2262 vs  4337
zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | grep -P '0[/|]1' | cut -f1,2,4,5 > $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.lst
#vcffilter -g "( GT = 0/1 )" $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | cut -f1,2,4,5 > $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.lst


M=$(wc -l $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.lst | cut -d" " -f1)
rm -f $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.info
for (( j = 1 ; j < $M+1 ; j++))
 do
	SCF="NC_014672.1"
	POS=$(cat $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.lst | sed -n $j'p' | cut -f2)
	#POS=198
	#zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | grep -n "0[/|]1"
	OCCURENCES=$(zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | grep -c "0/1")
	echo -n $SCF"\t"$POS"\t"$OCCURENCES"\t" >> $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.info
	for (( i = 1 ; i < $OCCURENCES+1 ; i++))
 		do
		SAMPLENR=$(zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -m 1 -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | grep -n "0/1" | sed -n $i'p' | cut -d : -f1)
		SAMPLENAME=$(vcfsamplenames $PROJECT/mt/vcf_master/mt_2.vcf.gz | sed -n $SAMPLENR'p')
		echo -n $SAMPLENAME"_" >> $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.info
		done
echo -n "\n" >> $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.info
echo $j"of"$M"\t"$SCF"\t"$POS
done

cat $PROJECT/mt/vcf_master/mt_2.heterozygous_positions.info | cut -f4 | tr '_' '\n' | sort | uniq -c | sort -nr > $PROJECT/mt/vcf_master/mt_2.heterozygous_individuals.info

cat $PROJECT/mt/vcf_master/mt_2.heterozygous_individuals.info
echo "SNPs $SNPS"
echo "SNPs without a heterozygous individual: $HMZ"
#3097 from 4350

############################################################################################################################################
##### this is not necessary when ploidy set to 1
######## htz to missing GT:GQ:DP:DPR:RO:QR:AO:QA:GL 0/1:-0:58:58,44:14:534:44:1705:-131.564,0,-30.8444 0/0:160.002:696:696,0:696:27261:0:0:0,-209.517,-2451.08
SCF="NC_014672.1"
POS=198
#POS=15548
# set Htz to missing
vcffilter -g "! ( GT = 0/1 )" $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100
#vice versa
vcffilter -g "( GT = 0/1 )" $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10| grep -n "0/1"
#vcffilter -s #filter entire records

vcffilter -g "(! GT = 0/1 )" $PROJECT/mt/vcf_master/mt_2.vcf.gz
vcffilter -g "( AO > RO )" $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100
zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | vcfbiallelic | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100
zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | vcfindels | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100
zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | vcfsnps | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100

#set missing genotype ( . ) to REF 0:.:.:.:.:.:.:.:.
vcffilter -g "! ( GT = 0/1 )" $PROJECT/mt/vcf_master/mt_2.vcf.gz | vcfglxgt -n | grep -v "#" | grep -m 1 -w -P $SCF'\t'$POS | tr '\t' '\n' | tail -n +10 | head -n 100

zcat $PROJECT/mt/vcf_master/mt_2.vcf.gz | grep -v "#" | cut -f 1,2,4,5,6
############################################################################################################################################

################# MT vcf
#mv $PROJECT/bin/MT.bed $PROJECT/bin/MT.temp
cat $PROJECT/bin/MT.temp |\
sed "s/MTinv_NC_014672/NC_014672.1/" |\
sed "s/__MTinv_NC_014672/__coding/" |\
sed "s/MT-B1__coding/MT-B1__PCRfragment/" |\
sed "s/MT-B2__coding/MT-B1__PCRfragment/" |\
sed "s/AT-rich__coding/AT-rich/" |\
sed "s/MT-ATrich__coding/MT-ATrich__PCRfragment/" > $PROJECT/bin/MT.bed

vcfannotate --bed $PROJECT/bin/MT.bed --key ANNOT --default NONCODING $PROJECT/mt/vcf_master/mt_2.vcf.gz |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.vcf.gz && \
tabix -fp vcf $PROJECT/mt/vcf_master/MT.vcf.gz

#################

######## filter indels
#zcat $PROJECT/mt/vcf_master/MT.vcf.gz | vcffilter -f "! ( TYPE = ins | TYPE = del | TYPE = complex )"
zcat $PROJECT/mt/vcf_master/MT.vcf.gz | vcffilter -f "! ( TYPE = ins | TYPE = del | TYPE = complex )" |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED | grep -v -P 'TYPE=ins|TYPE=del|TYPE=complex' |\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.noindel.vcf.gz
tabix -fp vcf $PROJECT/mt/vcf_master/MT.noindel.vcf.gz

#only biallelic
zcat $PROJECT/mt/vcf_master/MT.noindel.vcf.gz | vcfbiallelic | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED|\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.biallelic.vcf.gz
tabix -fp vcf $PROJECT/mt/vcf_master/MT.biallelic.vcf.gz

#only noindelsQ20
zcat $PROJECT/mt/vcf_master/MT.vcf.gz | vcffilter -f "QUAL > 20" | vcffilter -f "! ( TYPE = ins | TYPE = del | TYPE = complex )" | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
grep -v -P 'TYPE=ins|TYPE=del|TYPE=complex' |\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.noindelQ20.vcf.gz
tabix -fp vcf $PROJECT/mt/vcf_master/MT.noindelQ20.vcf.gz
#3822

#only biallelicQ20
zcat $PROJECT/mt/vcf_master/MT.noindel.vcf.gz | vcffilter -f "QUAL > 20" | vcfbiallelic | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED|\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz
tabix -fp vcf $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz
#3111

#zcat $PROJECT/mt/vcf_master/MT.biallelic.vcf.gz | grep -v "#" | cut -f1,2,4,5 | awk ' length($4) > 1'


### retain indels
zcat $PROJECT/mt/vcf_master/MT.vcf.gz | vcffilter -f "( TYPE = ins | TYPE = del | TYPE = complex )" |\
grep -v -P 'TYPE=snp|TYPE=mnp' | grep -v "TYPE=complex" |\
vcfuniq | bgzip -f -@ $CPUs -c /dev/stdin > $PROJECT/mt/vcf_master/MT.indels.vcf.gz
tabix -fp vcf $PROJECT/mt/vcf_master/MT.indels.vcf.gz


#Extract per sample
rm -rf $PROJECT/mt/vcf_individuals
mkdir -p $PROJECT/mt/vcf_individuals
vcfsamplenames $PROJECT/mt/vcf_master/MT.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.vcf.gz {} | vcffixup - | grep -v -P AO:QA:GL'\t'0: | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals/MT.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals/MT.{}.vcf.gz"
#4249

rm -rf $PROJECT/mt/vcf_individuals_noindels
mkdir -p $PROJECT/mt/vcf_individuals_noindels
vcfsamplenames $PROJECT/mt/vcf_master/MT.noindel.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.noindel.vcf.gz {} | vcffixup - | grep -v -P AO:QA:GL'\t'0: | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals_noindels/MT.noindel.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals_noindels/MT.noindel.{}.vcf.gz"
#4093

rm -rf $PROJECT/mt/vcf_individuals_biallelic
mkdir -p $PROJECT/mt/vcf_individuals_biallelic
vcfsamplenames $PROJECT/mt/vcf_master/MT.biallelic.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.biallelic.vcf.gz {} | vcffixup - | grep -v -P AO:QA:GL'\t'0: | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals_biallelic/MT.biallelic.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals_biallelic/MT.biallelic.{}.vcf.gz"
#3379

rm -rf $PROJECT/mt/vcf_individuals_noindelsQ20
mkdir -p $PROJECT/mt/vcf_individuals_noindelsQ20
vcfsamplenames $PROJECT/mt/vcf_master/MT.noindelQ20.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.noindelQ20.vcf.gz {} | vcffixup - | vcfbreakmulti | vcfuniq | grep -v -P AO:QA:GL'\t'0: | grep -v -P AO:QA:GL'\t''\.': | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals_noindelsQ20/MT.noindelQ20.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals_noindelsQ20/MT.noindelQ20.{}.vcf.gz"

#rm -rf $PROJECT/mt/vcf_individuals_biallelicQ20
mkdir -p $PROJECT/mt/vcf_individuals_biallelicQ20
vcfsamplenames $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz | parallel -j 1 "vcfkeepsamples $PROJECT/mt/vcf_master/MT.biallelicQ20.vcf.gz {} | vcffixup - | grep -v -P AO:QA:GL'\t'0: | bgzip -f -@ 30 -c /dev/stdin > $PROJECT/mt/vcf_individuals_biallelicQ20/MT.biallelicQ20.{}.vcf.gz; echo {}; tabix -fp vcf $PROJECT/mt/vcf_individuals_biallelicQ20/MT.biallelicQ20.{}.vcf.gz"





#################################################################


#################################################################


#################################################################


#################################################################



