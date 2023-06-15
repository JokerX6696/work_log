rm(list=ls())
setwd('D:\\desk\\XMSH_202305_4038')
library(grid)
library(ggplot2)
toBM <- read.table(file = 'Abundance_Stat.filter.anno.xls', sep = "\t", header = TRUE, encoding = 'UTF-8')


new_df <- toBM
# 属性矩阵
lab_df <- data.frame(toBM[,1:7])
names(lab_df) <- sub('\\.x','',names(lab_df))
# Cancer组 矩阵
cancer <- new_df[,grep('Tumor|HCC',names(new_df))]
# paracancer组 矩阵
paracancer <- new_df[,grep('Peritumor',names(new_df))]
for(num in seq(1,dim(cancer)[2]) ){
  names(cancer)[num] = paste('cancer','_',num,sep = '')
}
for(num in seq(1,dim(paracancer)[2]) ){
  names(paracancer)[num] = paste('paracancer','_',num,sep = '')
}
df <- data.frame(lab_df,cancer,paracancer)
pca_df <- df[,-c(1:7)]
rownames(pca_df) <- df$Specie

samples <- names(pca_df)
group <- c(rep('cancer',59),rep('paracancer',30))
mapping <- data.frame(Group=group)
rownames(mapping) <- samples


############################----------------------------------------------------------------

library(dplyr)
Phylum_df <- df[,-c(1,3,4,5,6,7)]

library(ggplot2)
ty = 'Genus'
cancer_2 <- apply(cancer,1,mean)
paracancer_2 <- apply(paracancer,1,mean)
heatmap_df <- data.frame(cancer=cancer_2,paracancer=paracancer_2)
#heatmap_df$mean <- apply(heatmap_df,1,mean)
heatmap_df$Genus <- df$Genus
df_sum <- heatmap_df %>%
  group_by(Genus) %>%
  summarize(across(c(cancer,paracancer),
                   sum, .names = "{.col}"))


#heatmap_df <- heatmap_df[rev(order(heatmap_df$mean)),]
#heatmap_df <- heatmap_df[1:30,-3]


# modified by tiansheng.xu; tiansheng.xu@oebiotech.com
library(ggplot2)
library("grid")
#library(oebio)

cols<-c("#FF0000FF","#FF9900FF", "#FFCC00FF" ,"#00FF00FF" ,"#6699FFFF", "#CC33FFFF",
        "#99991EFF" ,"#999999FF", "#FF00CCFF", "#CC0000FF", "#FFCCCCFF", "#FFFF00FF",
        "#CCFF00FF" ,"#358000FF", "#0000CCFF", "#99CCFFFF", "#00FFFFFF" ,"#CCFFFFFF",
        "#9900CCFF" ,"#CC99FFFF", "#996600FF", "#666600FF", "#666666FF" ,"#CCCCCCFF",
        "#79CC3DFF" ,"#CCCC99FF", "#63B8FF" ,  "#EECBAD"  , "#CD2626"   ,"#7CFC00","#99CCFFFF")

mycol<-cols
#mycol <-c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00","#CD2990","#9F79EE","#CD5C5C","#63B8FF","#EECBAD","#EEEE00","#7CFC00","#FF7F00","#C6E2FF")
###cols<-c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00","#CD2990","#9F79EE","#CD5C5C","#63B8FF","#EECBAD","#CD2626","#7CFC00","#FF7F00","#C6E2FF")
#a <- read.table(args[1],sep="\t",header=T,check.names=F)
###mycol<-c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00","#CD2990","#9F79EE","#CD5C5C","#63B8FF","#EECBAD","#CD2626","#7CFC00","#FF7F00","#C6E2FF")

otu1 <- as.data.frame(df_sum)#data.frame(Taxonomy=rownames(heatmap_df),heatmap_df)
counts <- apply(otu1[,2:ncol(otu1)],1,sum)
otu <- otu1[which(counts!=0),]
n<-ncol(otu)
#m<-log(n/2,2.85)
if(n<10){m=0.75
}else if(n>=10 & n<40){m<-log(n/2,9)
}else if(n>=40 & n<80){m<-log(n/2,7)
}else if(n>=80 & n<120){m<-log(n/2,3.5)
}else {m<-log(n/2,2.85)
}
if(nrow(otu)> 30){
  otu <- otu[1:30,]
}

if (nrow(otu)==1){q()}
if (nrow(otu)>length(mycol)) {
  library(RColorBrewer)
  cols<-c(cols,colorRampPalette(brewer.pal(8,"Dark2"))(nrow(otu)-length(cols)))
}
if (ncol(otu)>length(mycol)) {
  library(RColorBrewer)
  cols<-c(cols,colorRampPalette(brewer.pal(8,"Dark2"))(ncol(otu)-length(cols)))
}

rownames(otu) <- otu[,1]
otu <-otu[,-1]
al <- which(rownames(otu) %in% c("All"))
if(length(al)) otu <-otu[-al,]
rowsum <-sapply(1:nrow(otu),function(x) sum(otu[x,]))
otu<-otu[order(rowsum,decreasing=TRUE),]
dat <-sapply(1:ncol(otu),function(x) otu[,x]/sum(otu[,x]))
colnames(dat) <-colnames(otu)
rownames(dat) <-rownames(otu)
lab <-rownames(dat)
pdf(file=paste('Top30_',ty,'.pdf',sep = ''),height=8,width=9)
layout(matrix(1:2,2,1),heights=c(1.5:1))
par(mar=c(3,5,2,2))
barplot(dat*100,width=1,space=0.8,plot=T,las=1,col=mycol[1:nrow(dat)],cex.axis=1,cex.names=0.6/(m),border=NA,ylab="Relative abundance(%)",offset=0,cex.lab=1.5,las=2)
par(mar=c(2,2,1,0))
plot.new()
legend("topleft",legend=rownames(dat),ncol=2,fill=mycol[1:nrow(dat)],cex=0.8,bty="n",border=NA,x.intersp=0.2)
dev.off()














