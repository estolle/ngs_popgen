#!/usr/bin/env Rscript
#sumnumber.R $INPUTFILE
#sum numbers in a row

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
inputfile <- (args[1])

library(readr)
library(dplyr)

df <- read_tsv(inputfile, col_names = FALSE, skip = 1)
#drop.cols <- c(df$X1, df$X2)
#df %>% select_vars(-one_of(drop.cols)) -> df2
df %>% select (-c(X1, X2)) -> df2
df %>% select (c(X1, X2)) -> df3
Sum1 <- apply(df2, 1, function(x) sum(x==1,na.rm=TRUE))
Sum2 <- apply(df2, 1, function(x) sum(x==2,na.rm=TRUE))
Sum3 <- apply(df2, 1, function(x) sum(x==3,na.rm=TRUE))
data <- cbind(df3,Sum1,Sum2,Sum3)
#data2 <- sum(data$X3,data$X4,data$X5,na.rm=TRUE)
#data2 <- rowSums(data,[c(data$X3, data$X4, data$X5)], na.rm=TRUE)
data[,6] <- rowSums(data[,3:5])
write.table(data, file = paste0(inputfile,".zyg.sum.csv"), sep='\t', row.names=FALSE, col.names=FALSE)
