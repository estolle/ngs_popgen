#!/usr/bin/env Rscript
# E.Stolle MArch 2018

# supply with file (full path), KMER size and outputfolder
args <- commandArgs(trailingOnly = TRUE)
input <- (args[1])
KMER <- (args[2])
outputfolder <- (args[3])

library(findGSE)

#input <- "/scratch/ek/2018_Sol_MTreduce_kmer_repeat/depleted/40096151/kmc.35/CGI1Mir6-bigB-p.R1R2.trim.40096151.k35.kmc.count.txt.genomescopeformat.lst"

findGSE(histo=input, sizek=KMER, outdir=outputfolder)

print("done")

