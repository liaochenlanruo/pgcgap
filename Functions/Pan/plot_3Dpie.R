#!/usr/bin/env Rscript
# Take the output files from "--PAN" and plot pangenome 3D-pie.
# Author: Hualin Liu
#Date: April 8, 2019
#Version: 1.0.0

library("plotrix")

pdf("Pangenome_Pie.pdf")
stat = read.table("summary_statistics.txt", sep="\t")

type <- stat[1:4, 1]
def <- stat[1:4, 2]
num <- stat[1:4, 3]
total <- stat[5, 3]
per <- round(num/sum(num)*100,2)
percent <- paste(per, "%", sep="")

label <- paste(num,"(",percent,")", sep=" ")
labels <- paste(type, "\n", label)


i=1;
while(i<=4){
    if(num[i] <= 0){
        num[i] = 0.00001
    }else{
        num[i] = num[i]
    }
    i=i+1
}


par(oma=c(2,2,2,2))
pie3D(num, radius=0.75, height=0.15, labels=labels, explode=0.25, border="black", shade=0.6, cex.main=2, labelcex=1, main="Pangenome Pie")
dev.off()
