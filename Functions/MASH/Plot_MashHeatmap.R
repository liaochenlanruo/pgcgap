#!/usr/bin/env Rscript

require(pheatmap)
pdf("MASH_matrix.pdf")

mash <- read.delim("MASH.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)

if (nrow(mash) <= 20){
	wordsize = 15
}else if (20 < nrow(mash) & nrow(mash) <= 60){
	wordsize = 10
}else if (60 < nrow(mash) & nrow(mash) <= 80){
	wordsize = 8
}else if (80 < nrow(mash) & nrow(mash) <= 100){
	wordsize = 7
}else if (100 < nrow(mash) & nrow(mash) <= 145){
	wordsize = 5
}else if (145 < nrow(mash) & nrow(mash) <= 175){
	wordsize = 4
}else if (175 < nrow(mash) & nrow(mash) <= 250){
	wordsize = 3
}else if (250 < nrow(mash) & nrow(mash) <= 380){
	wordsize = 2
}else if (380 < nrow(mash) & nrow(mash) <= 700){
	wordsize = 1
}else if (700 < nrow(mash) & nrow(mash) <= 800){
	wordsize = 0.5
}else{
	wordsize = 0.1
}

aa <- pheatmap(mash, cluster_rows = TRUE, fontsize_row = wordsize, fontsize_col = wordsize)#added
order_row = aa$tree_row$order#get the row order
datat = data.frame(mash[order_row,order_row],check.names =F) #rearrange the original datas according to the new row order

if (any( mash == 0)){
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(100), colorRampPalette(colors = c("yellow","red"))(50)), legend_breaks = c(0, 75, 83, 90, 100), legend_labels = c("< 75", "75", "83", "90", "100"), cluster_rows = TRUE, cluster_cols = FALSE)
}else{
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, cluster_rows = TRUE, cluster_cols = FALSE)
}

dev.off()
