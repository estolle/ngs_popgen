#!/usr/bin/env Rscript
#QCplot.R $OUTPUTFOLDER trimming.stats

args <- commandArgs(trailingOnly = TRUE)
print(args[1])

work_dir <- (args[1])
inputfile <- (args[2])


setwd(work_dir)
#list.files(path = work_dir)

#library(extrafontdb)
library(extrafont)
#library(cowplot)
library(dplyr)
#library(tidyr)
library(ggplot2)
#library(tidyverse)
#library(RColorBrewer)
#library(scales)
#plyr
source('/home/xxxxx/scripts/es.r')
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


if(is.na(df2$DEGEN)) {
        df2$DEGEN <- 0
    }

if(is.na(df2$EMPTY)) {
        df2$EMPTY <- 0
    }

if(is.na(df2$QUALFILTERED)) {
        df2$QUALFILTERED <- 0
    }

if(is.na(df2$TOOSHORT)) {
        df2$TOOSHORT <- 0
    }


# Round to nearest even integer
#> round(-123.456,digits=-2)
#[1] -100
#ceiling(length(SAMPLES)) - ceiling(length(SAMPLES)) %% 2
#df$new <- 2 * round((length(SAMPLES))/2)

#even <- seq_len(length(SAMPLES)) %% 2   # index
#x.loadings <- data.frame(x=data.pc$loadings[even, ])
#y.loadings <- data.frame(x=data.pc$loadings[!even, ] )

#seq_len(length(SAMPLES)) %% 2
#factor(seq_len(length(SAMPLES)) %% 2)

data.frame(1:length(SAMPLES), grp1 = seq(from = 1, to = length(SAMPLES)), grp3=factor(seq_len(length(SAMPLES)) %% 2)) %>% mutate (grp2 = ifelse(grp3 == 0, "grey", "black")) %>% dplyr::filter(grp3 == 1) -> hlines
#data.frame(1:length(SAMPLES), grp1 = seq(from = 0.0, to = 0.001, by = 0.1), grp2 = c("grey","black"), grp3=factor(as.integer(df2$SAMPLE.f) %% 2)) %>% dplyr::filter(grp3 == 1) -> hlines
data.frame(1:length(SAMPLES), grp1 = seq(from = 1, to = length(SAMPLES)), grp3=factor(seq_len(length(SAMPLES)) %% 2)) %>% mutate (grp2 = ifelse(grp3 == 0, "grey", "black")) -> hlines_all
#data.frame(1:length(SAMPLES), grp1 = seq(from = 0.0, to = 0.001, by = 0.1),grp2= c("grey","black"), grp3 = factor(as.integer(df2$SAMPLE.f) %% 2)) -> hlines_all

title = ("reads_strict_afterTrim_afterDuplicates_raw")
ylabel <- title
MEANLINE <- df2.summary$READSIN_Mean
MEDIANLINE <- df2.summary$READSIN_Median
MEANLINE_after <- df2.summary$READPAIRSafter_Mean
MEDIANLINE_after <- df2.summary$READPAIRSafter_Mean
YMAX <- df2.summary$READSIN_Max
YMAX.f=as.integer(format((plyr::round_any(YMAX, 5000000, f = ceiling)), scientific=FALSE))


### points (with labels)
#scale_x_continuous(breaks=seq(from = 0, to = YMAX.f, by = (YMAX.f/5)), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f))
#annotate(geom="text", x = df2$READSIN, y = as.integer(df2$SAMPLE.f), label = ifelse(  ((df2$READSIN<(YMAX.f*0.4)) | (df2$READSIN>(YMAX.f*0.6))) ,as.character(sprintf("%.2f",df2$READSIN/1000000)),''), size = 3, colour = "black",hjust=-0.25, vjust=+0.5)
#geom_point(aes(x = df2$READPAIRSbefore, y = as.integer(df2$SAMPLE.f)), colour="black", fill="black", shape=21, size = 1.2, na.rm = FALSE)+
#geom_point(aes(x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+
#annotate(geom="text", x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f), label = paste(as.character(sprintf("%.1f",df2$SICKLEOUTPAIRED/1000000,'')),"__",df2$SICKLEPC,"%",sep=''), size = 2.5, colour = "darkgrey",hjust=-0.25, vjust=+0.5)
#geom_point(aes(x = df2$READPAIRSafter, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+

#scale_y_continuous




plot1 <-
    ggplot(df2, aes(READSIN, desc(SAMPLE.f),fill=factor(as.integer(SAMPLE.f) %% 2))) +
      geom_vline(xintercept = MEDIANLINE_after, linetype = 2, color="black",size=1.5) +
	geom_vline(xintercept = 0, linetype = 1, color="darkgrey",size=0.3) +
      geom_vline(xintercept = MEANLINE_after, linetype = 3, color="black",size=1.5) + 
      geom_vline(xintercept = MEDIANLINE, linetype = 2, color="darkgrey",size=0.5) +
      geom_vline(xintercept = MEANLINE, linetype = 3, color="darkgrey",size=0.5) + 
      geom_hline(data = hlines_all,aes(yintercept = hlines_all[,1]), size = 0.1, colour = "grey", linetype = 1 ,alpha=0.55) +
scale_x_continuous(breaks= c(1*(YMAX.f/10),1*(YMAX.f/3),2*(YMAX.f/3),3*(YMAX.f/3),MEDIANLINE_after), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f)) +
      scale_y_continuous(breaks = 1:length(SAMPLES), labels = function(y) {SAMPLES[y]},expand = c(0.1,0.1)) +
      scale_fill_manual(values = c('0' = "black", '1' = '#5599FF')) +
      labs(x="", y="", caption=title) +
      theme_es(grid='', legend.position='ylabel') +
      theme(axis.ticks.x = element_line(size=0.7), axis.ticks.y = element_line(size=0.7), axis.text.y=element_text(size = 10,margin = margin(0,6,0,0, unit = "pt")),axis.text.x=element_text(size = 10,margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt"),angle=-35)) +

geom_point(aes(x = 9*(YMAX.f/10), y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 3,na.rm = FALSE) +
geom_point(aes(x = 9*(YMAX.f/10)*1.02, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +
geom_point(aes(x = 9*(YMAX.f/10)*0.98, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.98, y = as.integer(length(df2$SAMPLE.f)), label = "untrimmed-remaining", size = 3, colour = "black",hjust=+1.5, vjust=+0.5,angle=90) +
annotate(geom="text", x = 9*(YMAX.f/10)*1.0, y = as.integer(length(df2$SAMPLE.f)), label = "after-filtering", size = 3, colour = "black",hjust=+1.5, vjust=+0.5,angle=90) +
annotate(geom="text", x = 9*(YMAX.f/10)*1.02, y = as.integer(length(df2$SAMPLE.f)), label = "before-filter", size = 3, colour = "black",hjust=+1.5, vjust=+0.5,angle=90) +

geom_point(aes(y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +
annotate(geom="text", x = df2$READSIN, y = as.integer(df2$SAMPLE.f), label = ifelse(  ((df2$READSIN<(YMAX.f*0.4)) | (df2$READSIN>(YMAX.f*0.6))) , as.character(sprintf("%.2f",df2$READSIN/1000000)),''), size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +
geom_point(aes(x = df2$READPAIRSafter, y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 3,na.rm = FALSE) +
annotate(geom="text", x = 1*(YMAX.f/10)*0.8, y = as.integer(df2$SAMPLE.f), label = as.character(sprintf("%.2f",df2$READPAIRSafter/1000000)), size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +
annotate(geom="text", x = 1*(YMAX.f/10)*1.5, y = as.integer(df2$SAMPLE.f), label = df2$READPAIRSafterPC, size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +
geom_point(aes(x = df2$READPAIRSafterUNTRIMMED, y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +


geom_point(x = df2$DEGEN, y = as.integer(df2$SAMPLE.f), colour="darkgrey", fill="darkgrey", shape=21, size = 2.0, na.rm = FALSE) +
geom_point(x = df2$QUALFILTERED/100, y = as.integer(df2$SAMPLE.f), colour="black", fill="black", shape=21, size = 1.0, na.rm = FALSE) +
geom_point(x = df2$TOOSHORT/100, y = as.integer(df2$SAMPLE.f), colour='black', fill='white', shape=21, size = 0.7, na.rm = FALSE) +
geom_point(x = df2$EMPTY/100, y = as.integer(df2$SAMPLE.f), colour="black", fill="black", shape=21, size = 0.3, na.rm = FALSE)


#loadfonts()
pdf(paste(inputfile,"_",title,".pdf",sep=""), width=29/2.54, height=21/2.54, useDingbats=FALSE)
print(plot1)
dev.off()



















