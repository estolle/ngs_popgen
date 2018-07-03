#!/bin/zsh
############################################################################
# reheader a vcf file (omitting scaffolds not scored)
# Eckart Stolle, Sept 2017
############################################################################
# usage: $HOME/scripts/vcfreheader.zsh INPUTVCF

if [ $# -ne 1 ]; then
    echo $0: usage: ./vcfreheader.zsh INPUTVCF
	echo "\nINPUTVCF: vcf file (bgzipped/tabix indexes)"
    exit 1
fi


INPUTVCF=$1
#INPUTVCF="$PROJECT/mt/vcf_master/mt_2.vcf.gz"
INPUTFILENAME=$(echo $INPUTVCF | rev | cut -d"." -f3- | rev)



tabix -H $INPUTVCF > $INPUTFILENAME.header.vcf

FIRSTLINE=$(echo $(grep -n "##contig" $INPUTFILENAME.header.vcf | head -n1 | cut -d":" -f1)-1 | bc -l) && echo $FIRSTLINE
LASTLINE=$(echo $(grep -n "##contig" $INPUTFILENAME.header.vcf | tail -n1 | cut -d":" -f1)+1 | bc -l) && echo $LASTLINE

zgrep -v "#" $INPUTVCF | cut -f1 | sort | uniq | sort -n | grep -v "\.1[A,B,C]$" > $INPUTFILENAME.header.scf_in_vcf.lst
zgrep "##contig=" $INPUTVCF > $INPUTFILENAME.header.scf_in_header.vcf

rm -f $INPUTFILENAME.header.new_scf.lst
N2=$(wc -l $INPUTFILENAME.header.scf_in_vcf.lst | cut -d' ' -f 1) && echo $N2
		for (( k = 1 ; k < $N2+1 ; k++))
		 do
	  		 CTG=$(cat $INPUTFILENAME.header.scf_in_vcf.lst | sed -n $k'p') && echo $CTG
	 	 	 grep $CTG $INPUTFILENAME.header.scf_in_header.vcf >> $INPUTFILENAME.header.new_scf.lst
	 	done

cat $INPUTFILENAME.header.vcf | head -n $FIRSTLINE > $INPUTFILENAME.newheader.vcf
cat $INPUTFILENAME.header.new_scf.lst >> $INPUTFILENAME.newheader.vcf
cat $INPUTFILENAME.header.vcf | tail -n +$LASTLINE >> $INPUTFILENAME.newheader.vcf

mv $INPUTVCF $INPUTVCF.bak
tabix -fp vcf -r $INPUTFILENAME.newheader.vcf $INPUTVCF.bak > $INPUTVCF
tabix -fp vcf $INPUTVCF

rm -f $INPUTFILENAME.header.new_scf.lst
rm -f $INPUTFILENAME.header.scf_in_header.vcf





