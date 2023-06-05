#Author: ycliu
#Date: 2017/03/23
rm(list = ls())
args<-commandArgs(T)

# setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
# args <-c("class.xls", "mapping.txt", "PCA_point_boxplot", "class.PERMANOVA.result.txt")
# args <-c("KEGG.level3.abundance.xls", "mapping.txt", "heatmap_KEGG")


if (length(args) != 4){
  cat('Usage: Rscript PCA_boxplot.r pyhlum.xls mapping.txt\n')
  q(status=1)
}
# print(args)
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

name <- gsub(".xls",'', args[1])
name <- gsub(".*/",'', name)
write.table(sum, paste0(args[3], "/",name,"_variance_of_component.txt"),quote=F,
            col.names = NA,
            row.names=T,sep="\t")

write.table(pc12_Group, paste0( args[3], "/",name,"_pc12_Group0.xls"),quote=F,sep="\t")

Data<-pca$x[,1:3]
#write.table(Data,args[4],quote=F,row.names=T,sep="\t")



fontsize1 = 17
fontsize2 = 19
linewidth = 1.2
Variance_PC1<-round(sum[2,1]*100,2)
Variance_PC2<-round(sum[2,2]*100,2)

theme_origin <- theme(text=element_text(family="serif"),  #设置字体
      panel.border = element_rect(size=linewidth, colour='black', fill=NA),
      plot.title = element_text(size = fontsize2,  #设置标题字体大小
                                colour = "black",  #设置标题颜色
                                # face = "bold",  #标题加粗
                                vjust = 0.7,
                                hjust = 0.5),  #设置标题居中
      
      axis.title.y = element_text(size = 17,  #设置y轴字体大小
                                  color = "black",  #设置y轴字体颜色
                                  #face = "bold",  #设置y轴标题加粗
                                  vjust = 1.9,  #调整y轴标题和轴线的距离
                                  hjust = 0.5,  #设置y轴标题居中
                                  angle = 90),  #设置y轴标题旋转角度
      
      axis.title.x = element_blank(),
      
      legend.title = element_text(color="black",  #修改图例的标题颜色
                                  size=17),  #图例的标题加粗
      
      legend.text = element_text(color="black",  #设置图例标签文字
                                 #face="bold",
                                 size = fontsize1),  #设置图例标签字体大小
      legend.position = c(1.05,0.7),
      # axis.text.x = element_blank(),
      axis.text.x = element_text(size = fontsize1,  #设置X轴标签字体大小
                                 color = "black",  #设置X轴标签颜色
                                 vjust = 1,  #设置X轴标签与轴线距离
                                 hjust = 0.5,  #设置X轴标签居中
                                 angle = 0),
      
      axis.text.y = element_text(size = fontsize1,  #设置y轴标签字体大小
                                 color = "black",  #设置y轴标签颜色
                                 vjust = 0.5,  #设置y轴标签与轴线距离
                                 hjust = 0.5,  #设置y轴标签居中
                                 angle = 0),  #设置y轴标签旋转角度
      
      axis.line.x = element_line(size=linewidth),  #设置X轴轴线大小
      axis.line.y = element_line(size=linewidth),  #设置X轴轴线大小
      axis.ticks.x= element_line(size=linewidth),  #设置X轴刻度线大小
      axis.ticks.y= element_line(size=linewidth))  ##设置y轴刻度线大小

PC1_min<- min(pc12_Group$PC1)
PC2_min<- min(pc12_Group$PC2)
PC1_max<- max(pc12_Group$PC1)
PC2_max<- max(pc12_Group$PC2)
PC1_scale <- (PC1_max-PC1_min)/100
PC2_scale <- (PC2_max-PC2_min)/100
PC1_x_lim = c(PC1_min-PC1_scale*20,PC1_max+PC1_scale*20)
#a <- ggplot(pc12_Group,aes(x=PC1,y=PC2,colour=Group))+geom_point(aes(shape=factor(Group)))+geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25)+xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+scale_x_continuous(limits=c(min,max+(max-min)/6))+theme_bw() +guides(col = guide_legend(nrow = 15))
a <- ggplot(pc12_Group,aes(PC1,PC2,group=Group,colour=Group))+
  geom_point(aes(PC1,PC2,colour=Group, shape=Group), size=3)+
  stat_ellipse(aes(x=PC1,y=PC2),size=linewidth,level = 0.6, type = "norm", show.legend = F)+
  # geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25,alpha=0.5)+
  # xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+
  scale_x_continuous(limits=PC1_x_lim, expand = c(0, 0))+
  # scale_y_continuous(limits=c(PC2_min-PC2_scale*5,PC2_max+PC2_scale*5), expand = c(0, 0))+
  ylab(paste("PC2 (",Variance_PC2,"%) ",sep=""))+
  theme_bw() +
  scale_colour_manual(values=mycol2[1:length(unique(pc12_Group$Group))]) + 
  guides(col = guide_legend(ncol = 1))+
  theme_origin+
  theme(
    plot.title = element_text(hjust = 0),
    axis.line.x = element_line(size=0, colour='black'),
    axis.line.y = element_line(size=0, colour='black'),
    panel.border = element_rect(size=linewidth, colour='black'),
    plot.margin = unit(c(0.6, 2, 2, 0.1), "inches"),
    panel.grid  =element_blank(),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    # legend.key.size = unit(5,"cm"),
    legend.key = element_blank(),
    legend.position = c(1.1,1.1),
    axis.text.x = element_text(color = NA),
    axis.ticks.x = element_line(color = NA)
  )

#P值
perP <- read.table(
  args[4],
  header = T, sep = "\t", check.names = F, row.names = 1,
  stringsAsFactors = FALSE
)
a <- a + labs(title=paste0("PERMANOVA Pvalue: ", perP$`Pr(>F)`[1]))

  # xlim(c(PC1_min-PC1_scale*5,PC1_max+PC1_scale*5))+
  # scale_x_reverse(limits=c(PC1_min-PC1_scale*5,PC1_max+PC1_scale*5), expand = c(0, 0))
#a <- ggplot(pc12_Group)+geom_point(aes(x=PC1,y=PC2,colour=Group),cex=2.5)+stat_ellipse(aes(x=PC1,y=PC2),level = 0.8, show.legend =F) +geom_text(aes(x=PC1,y=PC2,label=rownames(pc12_Group)),size=3,hjust=-0.1,vjust=1.25,alpha=0.5)+xlab(paste("PC1 ",Variance_PC1,"%",sep=""))+ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+scale_x_continuous(limits=c(min,max+(max-min)/6))+theme_bw() +guides(col = guide_legend(nrow = 15))
b <- ggplot(pc12_Group,aes(x=Group,y=PC1,fill=Group))+
  geom_boxplot(outlier.shape = NA,show.legend = F)+
  geom_jitter(aes(color=Group), shape = 21, size=3, alpha = 0.5, 
              colour='black',
              position = position_jitterdodge(jitter.width = 0.85,
                                              dodge.width = 0.75), #width =0.25,
              show.legend = F)+ ##添加散点
  scale_fill_manual(values=mycol2[1:length(unique(pc12_Group$Group))]) + 
  ylab(paste("PC1 (",Variance_PC1,"%) ",sep=""))+
  xlab('')+ #,axis.ticks.length = unit(0,"lines"),axis.ticks.margin = unit(1.65,"lines"))+xlab('')
  # ylim(c(PC1_min-PC1_scale*5,PC1_max+PC1_scale*5))+
  scale_y_continuous(limits=PC1_x_lim, expand = c(0, 0))+
  coord_flip()+
  # scale_y_reverse()+
  scale_x_discrete(breaks=levels(pc12_Group$Group),
                   labels=levels(pc12_Group$Group),
                   position='top')+
  theme_bw() +
  theme_origin+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.grid  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = fontsize2,  #设置X轴标签字体大小
                                       color = "black",  #设置X轴标签颜色
                                       vjust = 1,  #设置X轴标签与轴线距离
                                       hjust = 0.5,  #设置X轴标签居中
                                       angle = 0)
        )

p_pc1 <- round(t.test(PC1 ~ Group,data=pc12_Group)$p.value, 4)
b <- b+annotate("text", label=paste0("P = ", p_pc1), x=2.2,
                y= PC1_max+PC1_scale*-5,size = fontsize1*0.4,family='serif')

c <- ggplot(pc12_Group,aes(x=Group,y=PC2,fill=Group))+
  geom_boxplot(outlier.shape = NA,show.legend = F)+
  scale_fill_manual(values=mycol2[1:length(unique(pc12_Group$Group))]) + 
  geom_jitter(aes(color=Group), shape = 21, size=3, alpha = 0.5, 
              colour='black',
              position = position_jitterdodge(jitter.width = 0.85,
                                              dodge.width = 0.75), #width =0.25,
              show.legend = F)+ ##添加散点
  # scale_y_continuous(limits=c(PC2_min-PC2_scale*5,PC2_max+PC2_scale*5), expand = c(0, 0))+
  # ylab(paste("PC2 ",Variance_PC2,"%",sep=""))+
  # scale_y_continuous(limits=c(min,max+(max-min)/6))+
  # coord_flip()+
  theme_bw()+
  scale_x_discrete(breaks=levels(pc12_Group$Group),
                   labels=levels(pc12_Group$Group),
                   position='bottom')+
  theme(plot.margin = unit(c(0, 1, 0.5, 1.65), "lines"),
        axis.ticks.y = element_line(size = 0,color='white'))+
  xlab('')+#,axis.ticks.length = unit(0,"lines"),axis.ticks.margin = unit(1.65,"lines"))+xlab('')
  theme_origin+
  theme(
    plot.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid  = element_blank(),
    axis.line.y = element_blank(),
    # axis.line.y.left = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = fontsize1,  #设置X轴标签字体大小
                                color = "black",  #设置X轴标签颜色
                                vjust = 1,  #设置X轴标签与轴线距离
                                hjust = 0.5,  #设置X轴标签居中
                                angle = 0)
  )
p_pc2 <- round(t.test(PC2 ~ Group,data=pc12_Group)$p.value, 4)
c <- c+annotate("text", label=paste0("P = ", p_pc2), x=1.5,
                y= PC2_max-PC2_scale*5,size = fontsize1*0.4,family='serif')
# library(ggpubr)
# library(gridExtra)
# library(cowplot)
# p<-ggarrange(
#   ggarrange(a, c, widths=c(3:1), heights=c(1:1),align='h',
#             ncol=2, nrow=1), 
#   b, widths=c(3:1),heights=c(3:1),align='v',
#   ncol=1, nrow=2)


# p <- plot_grid(a, c, b, ncol = 2, nrow = 2,align = "hv",rel_heights = c(3,1),
#                rel_widths = c(3,2),
#                hjust=-1)
# 
# pdf(paste0(args[3], "/PCA_point_boxplot.pdf"), width = 20.4, height = 7.4)
# print(p)
# dev.off()

library(patchwork)
# p <- a+inset_element(c, left = 0.993, bottom = -0.183, right = 1.3, top = 1,align_to = "plot")
p <- a+inset_element(c, left = 0.992, bottom = 0, right = 1.3, top = 1,align_to = "plot")

p <- p+inset_element(b, left = 0, bottom = -0.5, right =1, top = 0,align_to = "panel")

# grid.newpage()
# pushViewport(viewport(layout=grid.layout(3,3),
#                       just="center"))
# vport<-function(x,y){
#   viewport(layout.pos.row=x,layout.pos.col=y)
# }
# 
pdf( paste0(args[3], "/",name,"_PCA_point_boxplot.pdf"))
print(p)
dev.off()

#write.table(sum,"variance_of_component.txt",quote=F,row.names=T,sep="\t")

