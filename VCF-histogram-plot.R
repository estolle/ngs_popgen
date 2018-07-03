#!/usr/bin/env Rscript
#Rscript scipt.R inputfile
#E. Stolle Jan 2018

args <- commandArgs(trailingOnly = TRUE)
print(args[1])

inputfile <- (args[1])

library(extrafont)
library(dplyr)
library(ggplot2)
library(patchwork)
source('/home/ek/scripts/es.r')

bold.16.text <- element_text(face = "bold", color = "black", size = 16)
bold.12.text <- element_text(face = "bold", color = "black", size = 12)
bold.10.text <- element_text(face = "bold", color = "black", size = 10)
bold.8.text <- element_text(face = "bold", color = "black", size = 8)
bold.6.text <- element_text(color = "black", size = 6)

df <- read_tsv(inputfile,col_names = FALSE)
#df <- read.csv(file=inputfile,header=FALSE, sep="\t",colClasses=c("character","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))
colnames(df) <- c("QUAL","AB","AC","AF","AN","AO","DP")

#df <- df %>% complete(SAMPLE, nesting(VALUE), fill = list(FRACTION = 0))
#SAMPLES <- levels(reorder(as.character(df$SAMPLE), desc(df$SAMPLE)))


df %>% summarise_if(.predicate = function(x) is.numeric(x), .funs = c(Mean="median", Max="max", Min="min"), na.rm = TRUE) -> df.summary

df %>% mutate(QUAL=replace(QUAL, QUAL>30000, 30000)) %>% mutate(AO=replace(AO, AO>10000, 10000)) %>% mutate(DP=replace(DP, DP>10000, 10000)) %>%
     as.data.frame() -> df2

QUALMAX <- df.summary$QUAL_Max
ABMAX <- df.summary$AB_Max
AFMAX <- df.summary$AF_Max
ANMAX <- df.summary$AN_Max
AOMAX <- df.summary$AO_Max
DPMAX <- df.summary$DP_Max
ACMAX <- df.summary$AC_Max


#works 
#+ geom_density(col="black")
#binwidth = 500
DP1 <- ggplot(df2, aes(df2$DP)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .1) +
  labs(x="DP (bin)", y="Count DP") +
  geom_freqpoly(bins = 200)
#  xlim(c(18,52)) + ylim(c(0,30))

DP2 <- ggplot(df2, aes(df2$DP)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .1) +
  labs(x="DP (bin)", y="Count DP") +
  geom_freqpoly(bins = 200) +
  xlim(c(0,3000))

AO1 <- ggplot(df2, aes(df2$AO)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .3) +
  labs(x="AO (bin)", y="Count AO") +
  geom_freqpoly(bins = 200)

AC1 <- ggplot(df2, aes(df2$AC)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .3) +
  labs(x="AC (bin)", y="Count AC") +
  geom_freqpoly(bins = 200)

AN1 <- ggplot(df2, aes(df2$AN)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .3) +
  labs(x="AN (bin)", y="Count AN") +
  geom_freqpoly(bins = 200)

AF1 <- ggplot(df2, aes(df2$AF)) +
  geom_histogram(bins = 100, col="black", fill="grey", alpha = .3) +
  labs(x="AF (bin)", y="Count AF") +
  geom_freqpoly(bins = 100)

AB1 <- ggplot(df2, aes(df2$AB)) +
  geom_histogram(bins = 100, col="black", fill="grey", alpha = .3) +
  labs(x="AB (bin)", y="Count AB") +
  geom_freqpoly(bins = 100)

QUAL1 <- ggplot(df2, aes(df2$QUAL)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .3) +
  labs(x="QUAL (bin)", y="Count QUAL") +
  geom_freqpoly(bins = 200)

QUAL2 <- ggplot(df2, aes(df2$QUAL)) +
  geom_histogram(bins = 200, col="black", fill="grey", alpha = .3) +
  labs(x="QUAL (bin)", y="Count QUAL") +
  geom_freqpoly(bins = 200) +
  xlim(c(0,10000))

combiplot <- QUAL1 / DP1 /
(QUAL2 | DP2) /
(AO1 | AN1 | AC1) /
(AB1 | AF1)

#loadfonts()
pdf(paste(inputfile,".vcf.stats.pdf",sep=""), width=21/2.54, height=29/2.54, useDingbats=FALSE)
print(combiplot)
dev.off()



















