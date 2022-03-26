#!/usr/bin/env Rscript

require(pheatmap)
pdf(file = "ANI_matrix.pdf", width = 8, height = 7)

ani <- read.delim("ANIs.heatmap", sep="\t", row.names=1, header=T, stringsAsFactors=FALSE, check.names=FALSE)
namecol <- names(ani)
maxlen <- unique(nchar(namecol[nchar(namecol)==max(nchar(namecol))]))

ani[is.na(ani)] <- 60
#wordsize = 4/nrow(ani)*50
cellwidth = 432/nrow(ani)
wordsize = 0.9*cellwidth/maxlen
if (wordsize < 2){
    wordsize = 2
}
miniani <- round(min(ani[ani>min(ani)]))
colorstepup <- 100 - miniani
colorsteplow <- miniani - 60

aa <- pheatmap(ani, cluster_rows = TRUE,fontsize_row = wordsize, fontsize_col = wordsize, silent = TRUE)#added
order_row = aa$tree_row$order#get the row order
datat = data.frame(ani[order_row,order_row],check.names =F) #rearrange the original datas according to the new row order
#datat[is.na(datat)] <- 60
if (any( ani == 60)){
	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(colorsteplow), colorRampPalette(colors = c("yellow","red"))(colorstepup)), legend_breaks = c(60, miniani, 95, 100), legend_labels = c("NA", miniani, "95", "100"), cluster_rows = TRUE, cluster_cols = FALSE, display_numbers = ifelse(datat > 95, "*", ""), cellwidth = cellwidth, cellhigh = cellwidth)
#	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(colorsteplow), colorRampPalette(colors = c("yellow","red"))(colorstepup)), legend_breaks = c(60, miniani, 85, 90, 95, 100), legend_labels = c("NA", miniani, "85", "90", "95", "100"), cluster_rows = TRUE, cluster_cols = FALSE, display_numbers = ifelse(datat > 95, "*", ""), cellwidth = cellwidth, cellhigh = cellwidth)
}else{
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, cluster_rows = TRUE, cluster_cols = FALSE, na_col = "#DDDDDD", display_numbers = ifelse(datat > 95, "*", ""), cellwidth = cellwidth, cellhigh = cellwidth)
}

if (any( ani == 60)){
	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(colorsteplow), colorRampPalette(colors = c("yellow","red"))(colorstepup)), legend_breaks = c(60, miniani, 95, 100), legend_labels = c("NA", miniani, "95", "100"), cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = cellwidth, cellhigh = cellwidth)
#	pheatmap(datat, border_color = "grey" ,fontsize_row = wordsize, fontsize_col = wordsize, color = c(colorRampPalette(colors = c("grey"))(colorsteplow), colorRampPalette(colors = c("yellow","red"))(colorstepup)), legend_breaks = c(60, miniani, 85, 90, 95, 100), legend_labels = c("NA", miniani, "85", "90", "95", "100"), cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = cellwidth, cellhigh = cellwidth)
}else{
	pheatmap(datat, border_color = "white" ,fontsize_row = wordsize, fontsize_col = wordsize, cluster_rows = TRUE, cluster_cols = FALSE, na_col = "#DDDDDD", cellwidth = cellwidth, cellhigh = cellwidth)
}
dev.off()
