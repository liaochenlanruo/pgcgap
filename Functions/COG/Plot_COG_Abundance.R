#!/usr/bin/env Rscript

library(gplots)
pdf("All_flags_relative_abundances.pdf")

cog <- read.table("All_flags_relative_abundances.table", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
cog2<-cog[,2:length(cog)]
cog3 <-data.matrix(cog2)
heatmap.2(cog3,col=bluered, labRow=row.names(cog),trace="none", scale="none", key=T,keysize=1.0,cexCol=1.2,cexRow=0.222, xlab="COG Functional Category", adjRow=c(0, 1))
dev.off()
