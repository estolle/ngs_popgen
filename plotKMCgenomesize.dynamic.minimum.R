#Rscript

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
inputfile <- (args[1])
TITLE <- (args[2])

kmercounts <- read.table(inputfile)
#plot(kmercounts[0:20,],type="l")
#points(kmercounts[2:20,])

#sum(as.numeric(kmercounts[2:10000,1]*kmercounts[2:10000,2]))
#2 103 603 758
#allKmers <- sum(as.numeric(kmercounts[2:10000,1]*kmercounts[2:10000,2]))
#max(kmercounts[2:60,2])
#45619475
peak <- max(kmercounts[3:1000,2])
#kmercounts$V1[kmercounts$V2==peak]
maxKmer <- kmercounts$V1[kmercounts$V2==peak]
peak.minimum <- min(kmercounts[0:(maxKmer),2])
minKmer <- kmercounts$V1[kmercounts$V2==peak.minimum]
allKmers <- sum(as.numeric(kmercounts[minKmer:10000,1]*kmercounts[minKmer:10000,2]))

ALLSUM <- round((allKmers/maxKmer)/1000000, digits = 2)
#ALLSUM
#525.9009 Mb

## only single copy region (peak between 2 and 15
#peakregionKmers <- sum(as.numeric(kmercounts[2:(10+maxKmer*3),1]*kmercounts[2:(10+maxKmer*3),2]))
peakregionKmers <- sum(as.numeric(kmercounts[minKmer:(10+maxKmer*3),1]*kmercounts[minKmer:(10+maxKmer*3),2]))
PEAKSUM <- round((peakregionKmers/maxKmer)/1000000, digits = 2)
#PEAKSUM
#357.3006 Mb
#--> n=1 kmers mostly seq error, but since peak is at n=4, its v close, thus the proximal end of the distribution is inside the n=1 peak, thus th$

# proportion
PROPORTION <- round(PEAKSUM/ALLSUM, digits = 4)
#PROPORTION
#0.6794067

REPETITIVE <- round(1-(PEAKSUM/ALLSUM), digits = 4)
#REPETITIVE
#0.3205933

PLOTTITLE1 <- paste(ALLSUM," Mb ",sep = "")
PLOTTITLE2 <- paste(PEAKSUM," Mb ",sep = "")
PLOTTITLE3 <- paste(REPETITIVE," repetitive",sep = "")

#For an inbred line, two set of genomes have mostly identical sequence and the major peak would be diploid (2C) peak. #On the otherhand, for high$

#compare the peak shape with poisson distribution
singleC <- (peakregionKmers/maxKmer)

png(paste(inputfile,".genomesize.png",sep=""),width = 1500, height = 1500)

plot(1:(20+maxKmer*3),dpois(1:(20+maxKmer*3), maxKmer)*singleC, type = "l", col=3, lty=2,xlab="kmer depth", ylab="kmers")
lines(kmercounts[1:(20+maxKmer*3),],type="l")
points(kmercounts[2:(20+maxKmer*3),])
points(kmercounts[maxKmer,],type = "p",pch = 16, col = "black",bg ="black", cex = 3)
points(kmercounts[minKmer,],type = "p",pch = 16, col = "black",bg ="black", cex = 3)
points(kmercounts[(10+maxKmer*3),],type = "p",pch = 16, col = "black",bg ="black", cex = 3)
text(maxKmer, (peak*1.01), labels=maxKmer, cex= 1,pos=3)
text(maxKmer+5, (peak*1.15), labels=paste(TITLE,PLOTTITLE1,PLOTTITLE2,PLOTTITLE3,sep = "\n"), cex= 2,pos=3)
dev.off()

new.df <- data.frame(NAME = c(TITLE), ALL = c(ALLSUM), PEAK = c(PEAKSUM), REP = c(REPETITIVE))
#write.table(new.df, file = paste0("/scratch/ek/testoutR",".summary.txt"), sep='\t', row.names=FALSE, col.names=FALSE)
write.table(new.df, file = paste0(inputfile,".summary.txt"), sep='\t', row.names=FALSE, col.names=FALSE)
print("done")

