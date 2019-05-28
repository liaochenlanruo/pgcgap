#!/usr/bin/env Rscript

require(pheatmap)
pdf("ANI_matrix.pdf")

ani <- read.delim("ANIs.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
pheatmap(ani, border_color = "white" ,fontsize_row = 3, fontsize_col = 3)

dev.off()
