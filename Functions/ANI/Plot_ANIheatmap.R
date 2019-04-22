#!/usr/bin/env Rscript

library(corrplot)
#tiff(file="Pangenome_matrix.tif", width=40, height=40, units="in", compression="lzw", res=150)
pdf("ANI_matrix.pdf")


ani <- read.delim("ANIs.heatmap", sep="\t", header=T, check.names=FALSE, stringsAsFactors=FALSE)
row.names(ani)<-ani$str
ani2 <- ani[,2:length(ani)]
ani3 <-data.matrix(ani2)

corrplot(ani3, method = "shade",shade.col = NA, tl.col ="black", tl.srt = 45, order = "AOE", is.corr = FALSE, cl.lim = c(80, 100), oma=c(2,2,2,2))

dev.off()
