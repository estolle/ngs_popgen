#!/usr/bin/env Rscript
#QCplot.R $OUTPUTFOLDER trimming.stats

args <- commandArgs(trailingOnly = TRUE)
print(args[1])

inputfile <- (args[1])

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
source('/home/ek/scripts/es.r')
#loadfonts()

bold.16.text <- element_text(face = "bold", color = "black", size = 16)
bold.12.text <- element_text(face = "bold", color = "black", size = 12)
bold.10.text <- element_text(face = "bold", color = "black", size = 10)
bold.8.text <- element_text(face = "bold", color = "black", size = 8)
bold.6.text <- element_text(color = "black", size = 6)



df <- read_tsv(inputfile,col_names = FALSE)
#df <- read.csv(file=inputfile,header=FALSE, sep="\t",colClasses=c("character","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))
colnames(df) <- c("SAMPLE","TOTALSNPs","MISSINGSNPsbefore","MISSINGSNPsafter","MISSINGSNPsincrease","TOTALheterzygous","TOTALhomozygous","InitiallyDetectedHeterozygous","TrueHeterozygous","FalseHeterozygousREF","FalseHeterozygousALT","multiallelicHeterozygousREFALT","multiallelicHeterozygousALTALT")

#df <- df %>% complete(SAMPLE, nesting(VALUE), fill = list(FRACTION = 0))
#SAMPLES <- levels(reorder(as.character(df$SAMPLE), desc(df$SAMPLE)))

df %>% mutate(SAMPLE.f = reorder(as.character(SAMPLE), desc(SAMPLE))) -> df2
df2 %>% summarise_if(.predicate = function(x) is.numeric(x), .funs = c(Mean="mean", Max="max",Median="median",Min="min"), na.rm = TRUE) -> df2.summary
SAMPLES <- levels(df2$SAMPLE.f)


#if(is.na(df2$DEGEN)) {
#        df2$DEGEN <- 0
#    }



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

title = ("VCF__filtering")
ylabel <- title
MEDIANMISSING <- df2.summary$MISSINGSNPsincrease_Median
MEDIANHTZ <- df2.summary$TOTALheterzygous_Median
MEDIANHMZ <- df2.summary$TOTALhomozygous_Median
YMAX <- df2.summary$TOTALSNPs_Max
YMAX.f=as.integer(format((plyr::round_any(YMAX, 1000000, f = ceiling)), scientific=FALSE))

Y1=as.integer(format((plyr::round_any(1*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y2=as.integer(format((plyr::round_any(2*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y3=as.integer(format((plyr::round_any(3*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y4=as.integer(format((plyr::round_any(4*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y5=as.integer(format((plyr::round_any(5*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y6=as.integer(format((plyr::round_any(6*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y7=as.integer(format((plyr::round_any(7*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y8=as.integer(format((plyr::round_any(8*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y9=as.integer(format((plyr::round_any(9*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))
Y10=as.integer(format((plyr::round_any(10*(YMAX.f/10), 1000000, f = ceiling)), scientific=FALSE))

### points (with labels)
#scale_x_continuous(breaks=seq(from = 0, to = YMAX.f, by = (YMAX.f/5)), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f))
#annotate(geom="text", x = df2$READSIN, y = as.integer(df2$SAMPLE.f), label = ifelse(  ((df2$READSIN<(YMAX.f*0.4)) | (df2$READSIN>(YMAX.f*0.6))) ,as.character(sprintf("%.2f",df2$READSIN/1000000)),''), size = 3, colour = "black",hjust=-0.25, vjust=+0.5)
#geom_point(aes(x = df2$READPAIRSbefore, y = as.integer(df2$SAMPLE.f)), colour="black", fill="black", shape=21, size = 1.2, na.rm = FALSE)+
#geom_point(aes(x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+
#annotate(geom="text", x = df2$SICKLEOUTPAIRED, y = as.integer(df2$SAMPLE.f), label = paste(as.character(sprintf("%.1f",df2$SICKLEOUTPAIRED/1000000,'')),"__",df2$SICKLEPC,"%",sep=''), size = 2.5, colour = "darkgrey",hjust=-0.25, vjust=+0.5)
#geom_point(aes(x = df2$READPAIRSafter, y = as.integer(df2$SAMPLE.f)), colour="darkgrey", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE)+





#scale_x_continuous(breaks= c(1*(YMAX.f/10),1*(YMAX.f/3),2*(YMAX.f/3),3*(YMAX.f/3),MEDIANLINE_after), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f)) +
#annotate(geom="text", x = df2$READSIN, y = as.integer(df2$SAMPLE.f), label = ifelse(  ((df2$READSIN<(YMAX.f*0.4)) | (df2$READSIN>(YMAX.f*0.6))) , as.character(sprintf("%.2f",df2$READSIN/1000000)),''), size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +
#annotate(geom="text", x = 1*(YMAX.f/10)*0.8, y = as.integer(df2$SAMPLE.f), label = as.character(sprintf("%.2f",df2$READPAIRSafter/1000000)), size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +
#annotate(geom="text", x = 1*(YMAX.f/10)*1.5, y = as.integer(df2$SAMPLE.f), label = df2$READPAIRSafterPC, size = 3, colour = "black",hjust=-0.25, vjust=+0.5) +

plot1 <-
ggplot(df2, aes(TOTALSNPs, desc(SAMPLE.f),fill=factor(as.integer(SAMPLE.f) %% 2))) +
      geom_vline(xintercept = 0, linetype = 1, color="darkgrey",size=0.3) +
      geom_vline(xintercept = MEDIANMISSING, linetype = 2, color="black",size=1.0) +
      geom_vline(xintercept = MEDIANHMZ, linetype = 2, color="black",size=1.0) + 
      geom_vline(xintercept = MEDIANHTZ, linetype = 2, color="black",size=1.0) +
      geom_hline(data = hlines_all,aes(yintercept = hlines_all[,1]), size = 0.1, colour = "grey", linetype = 1 ,alpha=0.55) +
    scale_x_continuous(breaks= c(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,MEDIANMISSING,MEDIANHTZ,MEDIANHMZ), labels = function(x) {sprintf("%00d", as.integer(x))}, expand = c(0,(YMAX.f/10)),limit = c(0, YMAX.f)) +
    scale_y_continuous(breaks = 1:length(SAMPLES), labels = function(y) {SAMPLES[y]},expand = c(0.1,0.1)) +
      scale_fill_manual(values = c('0' = "black", '1' = '#5599FF')) +
      labs(x="", y="", caption=title) +
      theme_es(grid='', legend.position='ylabel') +
      theme(axis.ticks.x = element_line(size=0.7), axis.ticks.y = element_line(size=0.7), axis.text.y=element_text(size = 8,margin = margin(0,6,0,0, unit = "pt")),axis.text.x=element_text(size = 10,margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt"),angle=-35)) +

geom_point(aes(x = 9*(YMAX.f/10)*1.07, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*1.07, y = as.integer(length(df2$SAMPLE.f)), label = "total", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +

geom_point(aes(x = 9*(YMAX.f/10)*1.05, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 2,na.rm = FALSE) +
geom_point(aes(x = 9*(YMAX.f/10)*1.04, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1.5,na.rm = FALSE) +
geom_point(aes(x = 9*(YMAX.f/10)*1.03, y = as.integer(length(df2$SAMPLE.f)), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*1.04, y = as.integer(length(df2$SAMPLE.f)), label = "missingPostPreDiff", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*1.00, y = as.integer(length(df2$SAMPLE.f)), colour="black", fill="darkgrey", shape=21 , size = 2.5) +
annotate(geom="text", x = 9*(YMAX.f/10)*1.00, y = as.integer(length(df2$SAMPLE.f)), label = "total HTZ", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*0.98, y = as.integer(length(df2$SAMPLE.f)), colour="white", fill="darkgrey", shape=21, size = 2.5,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.98, y = as.integer(length(df2$SAMPLE.f)), label = "total HMZ", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +

geom_point(x = 9*(YMAX.f/10)*0.95, y = as.integer(length(df2$SAMPLE.f)), colour="black", fill="darkgrey", shape=21, size = 1.0,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.95, y = as.integer(length(df2$SAMPLE.f)), label = "initial HTZ", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*0.93, y = as.integer(length(df2$SAMPLE.f)), colour="black", fill="darkgrey", shape=21, size = 1.0,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.93, y = as.integer(length(df2$SAMPLE.f)), label = "true HTZ", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*0.91, y = as.integer(length(df2$SAMPLE.f)), colour="black", fill="white", shape=21, size = 1.0,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.91, y = as.integer(length(df2$SAMPLE.f)), label = "false HTZ, REF", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*0.89, y = as.integer(length(df2$SAMPLE.f)), colour="black", fill="white", shape=21, size = 1.0,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.89, y = as.integer(length(df2$SAMPLE.f)), label = "false HTZ, ALT", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +

geom_point(x = 9*(YMAX.f/10)*0.86, y = as.integer(length(df2$SAMPLE.f)), colour='black', fill='black', shape=21, size = 0.3,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.86, y = as.integer(length(df2$SAMPLE.f)), label = "multiallelic-REF-ALT", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +
geom_point(x = 9*(YMAX.f/10)*0.85, y = as.integer(length(df2$SAMPLE.f)), colour='black', fill='black', shape=21, size = 0.3,na.rm = FALSE) +
annotate(geom="text", x = 9*(YMAX.f/10)*0.85, y = as.integer(length(df2$SAMPLE.f)), label = "multiallelic-ALT-ALT", size = 3, colour = "black",hjust=-0.2, vjust=-0.2,angle=-90) +

geom_point(aes(y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +

geom_point(aes(x = df2$MISSINGSNPsbefore, y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1.5,na.rm = FALSE) +
geom_point(aes(x = df2$MISSINGSNPsafter, y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 2,na.rm = FALSE) +
geom_point(aes(x = df2$MISSINGSNPsincrease, y = as.integer(SAMPLE.f), colour=factor(as.integer(SAMPLE.f) %% 2), fill=factor(as.integer(SAMPLE.f) %% 2)), size = 1,na.rm = FALSE) +

geom_point(x = df2$TOTALheterzygous, y = as.integer(df2$SAMPLE.f), colour="black", fill="darkgrey", shape=21, size = 2.5, na.rm = FALSE) +
geom_point(x = df2$TOTALhomozygous, y = as.integer(df2$SAMPLE.f), colour="white", fill="darkgrey", shape=21, size = 2.5, na.rm = FALSE) +

geom_point(x = df2$InitiallyDetectedHeterozygous, y = as.integer(df2$SAMPLE.f), colour="black", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE) +
geom_point(x = df2$TrueHeterozygous, y = as.integer(df2$SAMPLE.f), colour="black", fill="darkgrey", shape=21, size = 1.0, na.rm = FALSE) +
geom_point(x = df2$FalseHeterozygousREF, y = as.integer(df2$SAMPLE.f), colour="black", fill="white", shape=21, size = 1.0, na.rm = FALSE) +
geom_point(x = df2$FalseHeterozygousALT, y = as.integer(df2$SAMPLE.f), colour="black", fill="white", shape=21, size = 1.0, na.rm = FALSE) +

geom_point(x = df2$multiallelicHeterozygousREFALT, y = as.integer(df2$SAMPLE.f), colour='black', fill='black', shape=21, size = 0.3, na.rm = FALSE) +
geom_point(x = df2$multiallelicHeterozygousALTALT, y = as.integer(df2$SAMPLE.f), colour="black", fill="black", shape=21, size = 0.3, na.rm = FALSE)



#loadfonts()
pdf(paste(inputfile,".",title,".pdf",sep=""), width=21/2.54, height=29/2.54, useDingbats=FALSE)
print(plot1)
dev.off()



















