#!/usr/bin/env Rscript

library(ggplot2)

argv <- commandArgs(TRUE)
cog_tbl<- argv[1]
tbl<-read.table(cog_tbl,sep="\t",header=F)
names(tbl)<-c("id","flag","name")
p<-ggplot(data=tbl,aes(x=flag))

figure_name<-paste0(cog_tbl,".pdf")
pdf(figure_name,10,7)
par(oma=c(2,2,2,2))


p+geom_bar(aes(fill=flag),alpha=1)+
	theme_bw()+
	theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
#	scale_fill_manual(values=c(`A` = "#FCDCFC", `B` = "#FCDCCC", `C` = "#BCFCFC", `D` = "#FCFCDC", `E` = "#DCFCFC", `F` = "#DCECFC", `G` = "#CCFCFC", `H` = "#DCDCFC", `I` = "#DCCCFC", `J` = "#FCCCFC", `K` = "#FCDCEC", `L` = "#FCDCDC", `M` = "#ECFCAC", `N` = "#DCFCAC", `O` = "#9CFCAC", `P` = "#CCCCFC", `Q` = "#BCCCFC", `R` = "#E0E0E0", `S` = "#CCCCCC", `T` = "#FCFCAC", `U` = "#ACFCAC", `V` = "#FCFCBC", `W` = "#BCFCAC", `X` = "#9CFC9C", `Y` = "#FCFCCC", `Z` = "#CCFCAC"), guide = guide_legend(title = NULL, ncol=1), breaks=factor(unique(tbl$flag)), labels=unique(tbl$name))+
	scale_fill_manual(values=c(`A` = "#aec7e8", `B` = "#ff7f0e", `C` = "#ffbb78", `D` = "#2ca02c", `E` = "#98df8a", `F` = "#d62728", `G` = "#ff9896", `H` = "#9467bd", `I` = "#c5b0d5", `J` = "#8c564b", `K` = "#c49c94", `L` = "#e377c2", `M` = "#f7b6d2", `N` = "#c7c7c7", `O` = "#bcbd22", `P` = "#dbdb8d", `Q` = "#17becf", `R` = "#9edae5", `S` = "#393b79", `T` = "#9c9ede", `U` = "#b5cf6b", `V` = "#bd9e39", `W` = "#ccff00", `X` = "#ff33ff", `Y` = "#ffff00", `Z` = "#336666"), guide = guide_legend(title = NULL, ncol=1), breaks=factor(unique(tbl$flag)), labels=unique(tbl$name))+
	xlab("COG Functional Category")+
	ylab("Frequency")

invisible(dev.off())
