#!/usr/bin/env Rscript
# E.Stolle Febr 2018
# take from table regions according to bed-windows and for columns take means and some other calc
# Rscript medianwindows.R inputfile.data.csv inputfile.win.bed
# runs with 80% of all CPU cores!! takes relatively long

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
inputfile.data <- (args[1])
inputfile.win <- (args[2])

library(tidyverse)
library(parallel)

#library(foreach)
#library(doParallel)


#inputfile.data <- "/scratch/ek/apis/thelytoky/vcf/denise2017/popgen/geno.zyg.gz.zyg.sum.csv"
#inputfile.win <- "/scratch/ek/apis/thelytoky/vcf/denise2017/popgen/geno.50SNPs.25overlap.windows"

data <- read_tsv(inputfile.data, col_names = FALSE)
wins <- read_tsv(inputfile.win, col_names = FALSE)
#wins <- read_tsv(inputfile,na = "nan")

colnames(data) <- c("SCF","POS","Ns","HMZs","HTZs","SAMPLES")
colnames(wins) <- c("ID","START","END")

############################

wins.n <-as.integer(nrow(wins))
# 285426 lines/windows from the bed file

#add index to windows/bed file
wins$INDEX <- as.integer(c(1:(nrow(wins))))

#for(i in 1:wins.n)
#writeLines(c(SCFv,STARTv,ENDv,BASES,MED1,MED2,MED3,MED4,RAT1,RAT2,RAT3,RAT4),sep="\t")
#https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/

no_cores <- detectCores() * 0.65
#no_cores <- 105

cl <- makeCluster(no_cores, type="FORK")

system.time(
z <- do.call(rbind.data.frame,
parLapply(cl, 
          unique(wins$INDEX), 
          function(i){
	SCFv <- as.character(subset(wins, wins$INDEX == i, select=ID));
	STARTv <- as.numeric(subset(wins, wins$INDEX == i, select=START));
	ENDv <- as.numeric(subset(wins, wins$INDEX == i, select=END));
	BASES=ENDv-STARTv;
	subset.i <- subset(data, SCF == SCFv & POS >= STARTv & POS <= ENDv );
	#calculate Means
	MED1 <- round(mean(subset.i$Ns, na.rm = TRUE), digits = 4);
	MED2 <- round(mean(subset.i$HMZs, na.rm = TRUE), digits = 4);
	MED3 <- round(mean(subset.i$HTZs, na.rm = TRUE), digits = 4);
	MED4 <- round(mean(subset.i$SAMPLES, na.rm = TRUE), digits = 4);
	# ratio per sample
	RAT1 <- round(MED1/MED4, digits = 4);
	RAT2 <- round(MED2/MED4, digits = 4);
	RAT3 <- round(MED3/MED4, digits = 4);
	#ratio HTZ/HMZ
	RAT4 <- round(MED3/MED2, digits = 4);
      return(c(i,SCFv,STARTv,ENDv,BASES,MED1,MED2,MED3,MED4,RAT1,RAT2,RAT3,RAT4))}
))
)

stopCluster(cl)

#  user   system  elapsed 
#  45.127   34.732 4944.883
#1.5hrs


#1000SNPs
#user  system elapsed 
#0.149   0.120  21.340 

colnames(z) <- c("i","SCF","START","END","BASES","MED-N","MED-HMZ","MED-HTZ","MED-SAMPLES","RAT-NpSample","RAT-HMZpSample","RAT-HTZpSample","RAT-HTZpHMZ")

write.table(z, file = paste0(inputfile.win,".means.csv"), sep='\t', row.names=FALSE, col.names=TRUE)
print("done")

