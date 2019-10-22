#!/usr/bin/env Rscript
# Take the output files from "--PAN" and plot pangenome 3D-pie.
# Author: Hualin Liu
#Date: April 8, 2019
#Version: 1.0.0

library("plotrix")

pdf("Pangenome_Pie.pdf")
stat = read.table("summary_statistics.txt", sep="\t")

stats <- stat[rowSums(stat==0)==0,]
row <- nrow(stats)
row2 <- row-1

type <- stats[1:row2, 1]
def <- stats[1:row2, 2]
num <- stats[1:row2, 3]
total <- stats[row, 3]
per <- round(num/sum(num)*100,2)
percent <- paste(per, "%", sep="")

label <- paste(num,"(",percent,")", sep=" ")
labels <- paste(type, "\n", label)


par(oma=c(2,2,2,2))
pie3D(num, radius=0.75, height=0.1, labels=labels, explode=0, border="black", shade=0.6, cex.main=2, labelcex=0.9)
fan.plot(num, labels=labels)
dev.off()
