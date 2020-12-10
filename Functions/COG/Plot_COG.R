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

#p+geom_bar(aes(fill=flag),alpha=0.7)+
#	theme_bw()+
#	scale_fill_manual(values=rainbow(26), guide = guide_legend(title = NULL, ncol=1), breaks=factor(tbl$flag, levels=unique(as.character((tbl$flag)))), labels=tbl$name)+
#	xlab("Function Class")+
#	ylab("Frequency")

#mark <- unique(tbl$flag)
#lab <- unique(tbl$name)
#p+geom_bar(aes(fill=flag),alpha=0.7)+
#	theme_bw()+
#	scale_fill_manual(values=rainbow(26), guide = guide_legend(title = NULL, ncol=1), breaks=factor(mark), labels=lab)+
#	xlab("Function Class")+
#	ylab("Frequency")

p+geom_bar(aes(fill=flag),alpha=0.7)+
	theme_bw()+
	scale_fill_manual(values=rainbow(26), guide = guide_legend(title = NULL, ncol=1), breaks=factor(unique(tbl$flag)), labels=unique(tbl$name))+
	xlab("Function Class")+
	ylab("Frequency")

invisible(dev.off())
