#!/usr/bin/env Rscript

require(pheatmap)
pdf("ANI_matrix.pdf")

ani <- read.delim("ANIs.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(ani) <= 20){
	wordsize = 15
}else if (20 < nrow(ani) & nrow(ani) <= 60){
	wordsize = 10
}else if (60 < nrow(ani) & nrow(ani) <= 80){
	wordsize = 8
}else if (80 < nrow(ani) & nrow(ani) <= 100){
	wordsize = 7
}else if (100 < nrow(ani) & nrow(ani) <= 145){
	wordsize = 5
}else if (145 < nrow(ani) & nrow(ani) <= 175){
	wordsize = 4
}else if (175 < nrow(ani) & nrow(ani) <= 250){
	wordsize = 3
}else if (250 < nrow(ani) & nrow(ani) <= 380){
	wordsize = 2
}else if (380 < nrow(ani) & nrow(ani) <= 700){
	wordsize = 1
}else if (700 < nrow(ani) & nrow(ani) <= 800){
	wordsize = 0.5
}else{
	wordsize = 0.1
}
aa <- pheatmap(ani, cluster_rows = TRUE,fontsize_row = wordsize, fontsize_col = wordsize)#added
order_row = aa$tree_row$order#get the row order
datat = data.frame(ani[order_row,order_row],check.names =F) #rearrange the original datas according to the new row order

if (any( ani == 0)){
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(100), colorRampPalette(colors = c("yellow","red"))(50)), legend_breaks = c(0, 75, 83, 90, 100), legend_labels = c("< 75", "75", "83", "90", "100"), cluster_rows = TRUE, cluster_cols = FALSE)
}else{
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, cluster_rows = TRUE, cluster_cols = FALSE)
}

dev.off()
