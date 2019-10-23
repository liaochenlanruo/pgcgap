#!/usr/bin/env Rscript

require(pheatmap)
pdf("ANI_matrix.pdf")

ani <- read.delim("ANIs.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
aa <- pheatmap(ani, cluster_rows = TRUE)#added
order_row = aa$tree_row$order#get the row order
datat = data.frame(ani[order_row,order_row],check.names =F) #rearrange the original datas according to the new row order

if (any( ani == 0)){
	pheatmap(datat, border_color = "white" ,fontsize_row = 5, fontsize_col = 5, color = c(colorRampPalette(colors = c("grey"))(100), colorRampPalette(colors = c("yellow","red"))(50)), legend_breaks = c(0, 75, 83, 90, 100), legend_labels = c("< 75", "75", "83", "90", "100"), cluster_rows = TRUE, cluster_cols = FALSE)
}else{
	pheatmap(datat, border_color = "white" ,fontsize_row = 5, fontsize_col = 5, cluster_rows = TRUE, cluster_cols = FALSE)
}

dev.off()
