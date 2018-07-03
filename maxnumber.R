#!/usr/bin/env Rscript
#maxnumber.R $INPUTFILE
#max number in a row

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
inputfile <- (args[1])
data <- read.csv(file=inputfile, header=FALSE, sep=",",colClasses="integer")
colMax <- apply(data, 1, function(x) max(x,na.rm=TRUE))
data <- cbind(colMax, data)
write.csv(data, file = paste0(inputfile,".csv") ,row.names=FALSE)

