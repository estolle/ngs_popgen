#!/usr/bin/env Rscript
#QCplot.R $OUTPUTFOLDER trimming.stats

args <- commandArgs(trailingOnly = TRUE)
print(args[1])

work_dir <- (args[1])
inputfile <- (args[2])

#inputfile <- "trimming.stats"
#work_dir <- "/scratch/ek/reads/apis/Cridland2017/filtered_fastq/QC"
setwd(work_dir)
#list.files(path = work_dir)

library(extrafontdb)
library(extrafont)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(scales)
#plyr
source('/home/ek/scripts/es.r')
#loadfonts()

bold.16.text <- element_text(face = "bold", color = "black", size = 16)
bold.12.text <- element_text(face = "bold", color = "black", size = 12)
bold.10.text <- element_text(face = "bold", color = "black", size = 10)
bold.8.text <- element_text(face = "bold", color = "black", size = 8)
bold.6.text <- element_text(color = "black", size = 6)

df <- read_tsv(inputfile)
#colnames(df) <- c("SAMPLE","VALUE","FRACTION")
#df <- df %>% complete(SAMPLE, nesting(VALUE), fill = list(FRACTION = 0))
#SAMPLES <- levels(reorder(as.character(df$SAMPLE), desc(df$SAMPLE)))

df %>% mutate(SAMPLE.f = reorder(as.character(SAMPLE), desc(SAMPLE))) -> df2
df2 %>% summarise_if(.predicate = function(x) is.numeric(x), .funs = c(Mean="mean", Max="max",Median="median",Min="min"), na.rm = TRUE) -> df2.summary
SAMPLES <- levels(df2$SAMPLE.f)

data.frame(1:length(SAMPLES), grp1 = seq(from = 0.0, to = 0.001, by = 0.1),grp2= c("grey","black"), grp3 = factor(as.integer(df2$SAMPLE.f) %% 2)) %>% filter(grp3 == 1) -> hlines
data.frame(1:length(SAMPLES), grp1 = seq(from = 0.0, to = 0.001, by = 0.1),grp2= c("grey","black"), grp3 = factor(as.integer(df2$SAMPLE.f) %% 2)) -> hlines_all

title = ("reads_strict_afterTrim_afterDuplicates_raw")
ylabel <- title
MEANLINE <- df2.summary$READSIN_Mean
MEDIANLINE <- df2.summary$READSIN_Median
YMAX <- df2.summary$READSIN_Max
YMAX.f=as.integer(format((plyr::round_any(YMAX, 5000000, f = ceiling)), scientific=FALSE))


### points (with labels)
plot1 <-
    ggplot(df2, aes(READSIN, desc(SAMPLE.f),fill=factor(as.integer(SAMPLE.f) %% 2))) +
      geom_vline(xintercept = MEDIANLINE, linetype = 2) +
      geom_vline(xintercept = MEANLINE, linetype = 3) + 
      geom_hline(data = hlines_all,aes(yintercept = hlines_all[,1]), size = 0.1, colour = "grey", linetype = 1 ,alpha=0.55) +
      scale_x_continuous(breaks=seq(from = 0, to = YMAX.f, by = (YMAX.f/5)), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f)) +
      scale_y_continuous(breaks = 1:length(SAMPLES), labels = function(y) {SAMPLES[y]},expand = c(0.1,0.1)) +
      scale_fill_manual(values = c('0' = "black", '1' = '#5599FF')) +
      labs(x="", y="", caption=title) +
      theme_es(grid='', legend.position='ylabel') +
      theme(axis.ticks.x = element_line(size=0.7), axis.ticks.y = element_line(size=0.7), axis.text.y=element_text(size = 10,margin = margin(0,6,0,0, unit = "pt")),axis.text.x=element_text(size = 10,margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt"))) +

geom_point(aes(y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 3,na.rm = FALSE) +
annotate(geom="text", x = df2$READSIN, y = as.integer(df2$SAMPLE.f), label = ifelse(  ((df2$READSIN<(YMAX.f*0.4)) | (df2$READSIN>(YMAX.f*0.6))) ,as.character(sprintf("%.2f",df2$READSIN/1000000)),''), size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +

geom_point(aes(x = df2$READPAIRSbefore, y = as.integer(df2$SAMPLE.f)), colour="black", fill="black", shape=21, size = 1.2, na.rm = FALSE)+
geom_point(aes(x = df2$READPAIRSafter, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+
geom_point(aes(x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+
annotate(geom="text", x = 0, y = as.integer(df2$SAMPLE.f), label = df2$READPAIRSafterPC, size = 3, colour = "black",hjust=+1.5, vjust=+0.5)+
annotate(geom="text", x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f), label = paste(as.character(sprintf("%.1f",df2$SICKLEOUTPAIRED/1000000,'')),"__",df2$SICKLEPC,"%",sep=''), size = 2.5, colour = "darkgrey",hjust=-0.25, vjust=+0.5)

#loadfonts()
pdf(paste(inputfile,"_",title,".pdf",sep=""), width=29/2.54, height=21/2.54, useDingbats=FALSE)
print(plot1)
dev.off()

###########################
title = ("filtered")
ylabel <- title
#df2$SICKLEOUTSINGLE/1000
#df2$OPTDUPL
#df2$DEGEN
#df2$QUALFILTERED/100
#df2$TOOSHORT/100
#df2$EMPTY

YMAX=60000
YMAX.f=YMAX

plot2 <-
    ggplot(df2, aes((SICKLEOUTSINGLE/1000), desc(SAMPLE.f),fill=factor(as.integer(SAMPLE.f) %% 2))) +
      geom_hline(data = hlines_all,aes(yintercept = hlines_all[,1]), size = 0.1, colour = "grey", linetype = 1 ,alpha=0.55) +
      scale_x_continuous(breaks=seq(from = 0, to = YMAX.f, by = (YMAX.f/5)), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f)) +
      scale_y_continuous(breaks = 1:length(SAMPLES), labels = function(y) {SAMPLES[y]},expand = c(0.1,0.1)) +
      scale_fill_manual(values = c('0' = "black", '1' = '#5599FF')) +
      labs(x="", y="", caption=title) +
      theme_es(grid='', legend.position='ylabel') +
      theme(axis.ticks.x = element_line(size=0.7), axis.ticks.y = element_line(size=0.7), axis.text.y=element_text(size = 10,margin = margin(0,6,0,0, unit = "pt")),axis.text.x=element_text(size = 10,margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt"))) +

geom_point(aes(y = as.integer(SAMPLE.f), colour='#5599FF', fill='#5599FF'), size = 3,na.rm = FALSE) +
geom_point(aes(x = df2$OPTDUPL, y = as.integer(df2$SAMPLE.f)), colour="black", fill="black", shape=21, size = 1.2, na.rm = FALSE)+
geom_point(aes(x = df2$DEGEN, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+
geom_point(aes(x = df2$QUALFILTERED/100, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=20, size = 1.0, na.rm = FALSE)+
geom_point(aes(x = df2$TOOSHORT/100, y = as.integer(df2$SAMPLE.f)), colour='#5599FF', fill="black", shape=19, size = 1.0, na.rm = FALSE)+
geom_point(aes(x = df2$EMPTY/100, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="black", shape=18, size = 1.0, na.rm = FALSE)+

geom_point(aes(x = YMAX*0.9, y = length(SAMPLES), colour='#5599FF'), fill='#5599FF', size = 3) +
geom_point(aes(x = YMAX*0.9, y = length(SAMPLES)-1, colour="black"), fill="black", shape=21, size = 1.2)+
geom_point(aes(x = YMAX*0.9, y = length(SAMPLES)-2, colour="darkgrey"), fill="darkgrey", shape=21, size = 1.0)+
geom_point(aes(x = YMAX*0.9, y = length(SAMPLES)-3, colour="darkgrey"), fill="darkgrey", shape=20, size = 1.0)+
geom_point(aes(x = YMAX*0.9, y = length(SAMPLES)-4, colour='#5599FF'), fill="black", shape=19, size = 1.0)+
geom_point(aes(x = YMAX*0.9, y = length(SAMPLES)-5, colour="darkgrey"), fill="black", shape=18, size = 1.0)+

annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES), size = 3, hjust = 0, label = "singletons(reads)(highQ)(*1000)")+
annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES)-1, size = 3, hjust = 0, label = "optical duplicates")+
annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES)-2, size = 3, hjust = 0, label = "N-containing reads")+
annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES)-3, size = 3, hjust = 0, label = "low quality reads(*100)")+
annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES)-4, size = 3, hjust = 0, label = "too short reads(*100)")+
annotate(geom = "text", x = YMAX*0.92, y = length(SAMPLES)-5, size = 3, hjust = 0, label = "empty readpairs")

#loadfonts()
pdf(paste(inputfile,"_",title,".pdf",sep=""), width=29/2.54, height=21/2.54, useDingbats=FALSE)
print(plot2)
dev.off()

