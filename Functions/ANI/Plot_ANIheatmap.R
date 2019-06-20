#!/usr/bin/env Rscript

require(pheatmap)
pdf("ANI_matrix.pdf")

ani <- read.delim("ANIs.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
if (any( ani == 0)){
	pheatmap(ani, border_color = "white" ,fontsize_row = 3, fontsize_col = 3, color = c(colorRampPalette(colors = c("grey"))(100), colorRampPalette(colors = c("yellow","red"))(50)), legend_breaks = c(0, 75, 83, 90, 100), legend_labels = c("< 75", "75", "83", "90", "100"))
}else{
	pheatmap(ani, border_color = "white" ,fontsize_row = 3, fontsize_col = 3)
}

dev.off()
