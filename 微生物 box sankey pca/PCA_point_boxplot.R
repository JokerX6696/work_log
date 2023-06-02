#Author: ycliu
#Date: 2017/03/23
rm(list = ls())
args<-commandArgs(T)

setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
args <-c("class.xls", "mapping.txt", "PCA_point_boxplot")
# args <-c("KEGG.level3.abundance.xls", "mapping.txt", "heatmap_KEGG")


if (length(args) != 3){
  cat('Usage: Rscript PCA_boxplot.r pyhlum.xls mapping.txt\n')
  q(status=1)
}

library(ggplot2)
library(grid)
# library("plot3D")
# library("scatterplot3d")
# library("maptools")


abundence <-read.table(args[1],header=T,sep="\t",check.names=F,row.names=1)
mapping  <-read.table(args[2],header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=FALSE,comment.char="")#modified by xtsheng 20180424
fa<- unique(mapping$Group)#modified by xtsheng 20180424
mapping$Group<-factor(mapping$Group,levels=fa)#modified by xtsheng 20180424
n<-length(colnames(abundence))

mycol2<-rep(c( "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#F0A3FF","#740AFF","#990000","#FFFF00","#CD2990","#9F79EE","#CD5C5C","#63B8FF","#EECBAD","#CD2626","#7CFC00","#FF7F00","#C6E2FF"),n)
mypch<-rep(c("16","17","18","19","15","9","4","3","23","22","25","7"),n)
mycol2 <- c("#0072B5", "#BC3C29")


abundence<-abundence[apply(abundence,1,var) != 0,]
pca <- prcomp(t(abundence),cor=T,sca=T)

pc12<-pca$x[,1:2]
pc12<-as.data.frame(pc12)
pc12_Group<-merge(pc12,mapping[,'Group',drop=F],by.x='row.names',by.y='row.names')

max<-round(max(pc12))
min<-min(pc12)

rownames(pc12_Group)<-pc12_Group[,1]
pc12_Group<-pc12_Group[,-1]

sum<-summary(pca)
sum<-as.data.frame(sum$importance)

if(!dir.exists(args[3])){
  dir.create(args[3])
}

write.table(sum, paste0(args[3], "/variance_of_component.txt"),quote=F,
            col.names = NA,
            row.names=T,sep="\t")

write.table(pc12_Group, paste0( args[3], "/pc12_Group0.xls"),quote=F,sep="\t")

Data<-pca$x[,1:3]
#write.table(Data,args[4],quote=F,row.names=T,sep="\t")




Variance_PC1<-round(sum[2,1]*100,2)
Variance_PC2<-round(sum[2,2]*100,2)


#a <- ggplot(pc12_Group,aes(x=PC1,y=PC2,colour=Group))+geom_point(aes(shape=factor(Group)))+geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25)+xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+scale_x_continuous(limits=c(min,max+(max-min)/6))+theme_bw() +guides(col = guide_legend(nrow = 15))
a <- ggplot(pc12_Group,aes(PC1,PC2,group=Group,colour=Group))+geom_point(aes(PC1,PC2,colour=Group))+stat_ellipse(aes(x=PC1,y=PC2),level = 0.6, show.legend =F)+geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25,alpha=0.5)+xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+theme_bw() +guides(col = guide_legend(nrow = 15))
#a <- ggplot(pc12_Group)+geom_point(aes(x=PC1,y=PC2,colour=Group),cex=2.5)+stat_ellipse(aes(x=PC1,y=PC2),level = 0.8, show.legend =F) +geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25,alpha=0.5)+xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+scale_x_continuous(limits=c(min,max+(max-min)/6))+theme_bw() +guides(col = guide_legend(nrow = 15))
b <- ggplot(pc12_Group,aes(x=Group,y=PC1,fill=Group))+geom_boxplot()+ylab(paste("PC1 ",Variance_PC1,"%",sep=""))+scale_y_continuous(limits=c(min,max+(max-min)/6))+coord_flip()+theme_bw()+scale_x_discrete(breaks=pc12_Group$Group,labels=rep('',length(pc12_Group$Group)))+theme(plot.margin = unit(c(0, 1, 0.5, 1.65), "lines"),axis.ticks.y = element_line(size = 0,color='white'))+xlab('')#,axis.ticks.length = unit(0,"lines"),axis.ticks.margin = unit(1.65,"lines"))+xlab('')
c <- ggplot(pc12_Group,aes(x=Group,y=PC2,fill=Group))+geom_boxplot()+ylab(paste("PC1 ",Variance_PC1,"%",sep=""))+scale_y_continuous(limits=c(min,max+(max-min)/6))+coord_flip()+theme_bw()+scale_x_discrete(breaks=pc12_Group$Group,labels=rep('',length(pc12_Group$Group)))+theme(plot.margin = unit(c(0, 1, 0.5, 1.65), "lines"),axis.ticks.y = element_line(size = 0,color='white'))+xlab('')#,axis.ticks.length = unit(0,"lines"),axis.ticks.margin = unit(1.65,"lines"))+xlab('')

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))
vport<-function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)
}

pdf( paste0(args[3], "/PCA_circo_boxplot.pdf"))
print(a,vp=vport(1:2,1))
print(b,vp=vport(3,1))
dev.off()

#write.table(sum,"variance_of_component.txt",quote=F,row.names=T,sep="\t")

