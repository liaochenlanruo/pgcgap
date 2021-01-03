#!/usr/bin/env Rscript

require(pheatmap)
pdf(file = "MASH_matrix.pdf", width = 8, height = 7)

mash <- read.delim("MASH.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)

wordsize = 7/nrow(mash)*50

aa <- pheatmap(mash, cluster_rows = TRUE, fontsize_row = wordsize, fontsize_col = wordsize, silent = TRUE)#added
order_row = aa$tree_row$order#get the row order
datat = data.frame(mash[order_row,order_row],check.names =F) #rearrange the original datas according to the new row order

#if (any( mash == 0)){
#	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(100), colorRampPalette(colors = c("yellow","red"))(50)), legend_breaks = c(0, 75, 83, 90, 100), legend_labels = c("< 75", "75", "83", "90", "100"), cluster_rows = TRUE, cluster_cols = FALSE)
#}else{
	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, cluster_rows = TRUE, cluster_cols = FALSE, display_numbers = ifelse(datat > 95, "*", ""))
#}

dev.off()
