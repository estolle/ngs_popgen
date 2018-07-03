mkdir -p ~/test
#-lib $TEs -dir ~/test
REF="/scratch/genomes/genomes_original/Amel/chromosome_assembly_GCA_000002195.1/GCA_000002195.1_Amel_4.5_genomic.fna"
RepeatMasker -div 10 -gff -u -xm -html -xsmall -lcambig -a -frag 60000 -nolow -qq -pa 110 -cutoff 350 -dir ~/test $REF

