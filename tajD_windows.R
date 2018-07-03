#!/usr/bin/env Rscript
# E.Stolle Febr 2018
# Rscript tajD_windows.R inputfile.data.tsv


args <- commandArgs(trailingOnly = TRUE)
print(args[1])
inputfile.data <- (args[1])


library(tidyverse)
library(parallel)
library(dplyr)


df <- read_tsv(inputfile.data, col_names = TRUE)
#na = "nan"
#colnames(data) <- c("SCF","POS","Ns","HMZs","HTZs","SAMPLES")
#wins.n <-as.integer(nrow(df))
# 285426 lines/windows

#add index to windows/bed file
df$INDEX <- as.integer(c(1:(nrow(df))))

######## define TajD function

#a simple function to calculate Tajima's D from nucleotide diversity pi, the number of segregating sites, and the number of sequences/samples
#https://github.com/thomasblankers/popgen/blob/master/TajimaD.r

harmonic<-function(x) { out<-c(); for(i in 1:(x-1)) { out<-c(out,1/i)}; sum_out=sum(out); sum_out} #calculate the (n-1)th harmonic number
harmonic_sq<-function(x) { out<-c(); for(i in 1:(x-1)) { out<-c(out,1/(i^2))}; sum_out=sum(out); sum_out}
TajimaD <- function(pi,S,n) {
	a1=harmonic(n)
	a2=harmonic_sq(n)
	b1=(n+1)/(3*(n-1))
	b2=(2*((n^2)+n+3))/((9*n)*(n-1))
	c1=b1-(1/a1)
	c2=b2-((n+2)/(a1*n))+(a2/(a1^2))
	e1=c1/a1
	e2=c2/((a1^2)+a2)
	d=(pi-(S/a1)) # calculate difference between observed nucleotide diversity and expected nucleotide diversity
	D=d/sqrt(e1*S+(e2*S)*(S-1)) # calculate Tajima's D
	D
}


#####################################
## Calculate TajD from cols 6(pi1),5(SegrSites),14(Samples)

#data2$taj <- apply(data2, 1, function(x) TajimaD(x[1],x[2],x[3]))
### works #42 sec for 280k rows

system.time(
df$taj1 <- apply(df[,c(6,5,14)], 1, function(x) TajimaD(x[1],x[2],x[3]))
)
system.time(
df$taj2 <- apply(df[,c(7,5,14)], 1, function(x) TajimaD(x[1],x[2],x[3]))
)


#####################################
# ratio between 2 different pi values

df <- mutate(df, pi_diff = pi_th - pi_ar)

#####################################
# ratio between 2 different TajD values

df <- mutate(df, taj_diff = taj2 - taj1)

write.table(df, file = paste0(inputfile.data,".with_tajD.tsv"), sep='\t', row.names=FALSE, col.names=TRUE)
print("done")

