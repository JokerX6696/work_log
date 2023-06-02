rm(list=ls())

# 并行初始化
library(ggplot2)
library(reshape2)
library(ggimage)
args<-commandArgs(T)

# setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
# args<- c("L2.xls",
#          "mapping.txt",
#          "diff_boxplot",
#          "t")

if (length(args) != 4){
  cat('Usage: anova_post.hoc.r data.txt mapping.txt results.txt heatmap.txt aov\n aov: anova\n kw: kruskal wallis\n t: t test\n w: wilcoxon \n')
  q(status=1)
}


filedata <- read.table(args[1],sep="\t",header=T,row.names=1,check.names=F, quote = "")
mapping <- read.table(args[2], sep="\t",header=F,check.names=F)
colnames(mapping) <- mapping[1,]
mapping <- mapping[-1,]
row.names(mapping) <- mapping$Sample
num <- length(rownames(filedata))


#1.计算P值

res <- NULL
for(i in 1:num){
  if(sd((filedata[i,]))!=0){
    x <- as.data.frame(t(filedata[i,]))
    colnames(x) <- c("datainfo")
    merged_data<-merge(x,mapping[,'Group',drop=F],by.x='row.names',by.y='row.names')
    if(args[4]=="aov"){
      method <- "ANOVA"
      merged_data.aov <- aov(datainfo ~ Group,data=merged_data)
      pvalue <- summary(merged_data.aov)[[1]][5][1,]
    }else if(args[4]=="kw"){
      method <- "Kruskal Wallis"
      p <- kruskal.test(datainfo ~ Group,data=merged_data)
      pvalue <- p$p.value
    }else if(args[4]=="t"){
      method <- "T test"
      
      t <- 0
      for (j in 1:length(unique(as.character(merged_data$Group)))) {
        a <- merged_data[merged_data$Group == unique(as.character(merged_data$Group))[j],]
        if (length(unique(a$datainfo)) == 1) t <- t + 1
      }
      
      if (t == length(unique(as.character(merged_data$Group)))) {
        message("组内存在恒量: ", x)
        next
      }
      
      p <- t.test(datainfo ~ Group,data=merged_data)
      pvalue <- p$p.value
    }else if(args[4]=="w"){
      method <- "Wilcoxon"
      p <- wilcox.test(datainfo ~ Group,data=merged_data)
      pvalue <- p$p.value
    }
    
    if(pvalue<=0.001){
      significant <- "***"
    }else if(pvalue > 0.001 && pvalue<= 0.01){
      significant <- "**"
    }else if(pvalue > 0.01 && pvalue<= 0.05){
      significant <- "*"
    }else if(pvalue>0.05){
      significant <- "not significant"
    }
    
    df <- cbind(filedata[i,],pvalue,significant)
    res <- rbind(res,df)
  }
}


#####
Index <- rownames(res)
results <- cbind(Index,as.data.frame(res))
results <- as.data.frame(results)
results$FDR_P <- p.adjust(as.numeric(as.character(results$pvalue)), "BH")
# heatmap <- results[as.numeric(as.character(results$pvalue))<0.05,]
# heatmap <- subset(heatmap, select = -c(pvalue, significant, FDR_P))
name <- gsub("^.*L", "", args[1])
name <- gsub(".xls", "", name)

if(!dir.exists(args[3])){
  dir.create(args[3])
}
write.table(results,paste0(args[3], "/L",name ,".diff.xls"),sep="\t",quote=F,row.names=F)

#取top30

results1 <- cbind(results, do.call(rbind, strsplit(row.names(results), ";")))
results1[,name] <- do.call(rbind, strsplit(results1[,name], "__"))[, 2]
results1[,'1'] <- do.call(rbind, strsplit(results1[,'1'], "__"))[, 2]
results1 <- results1[!grepl("unclassified", results1[, name]), ] #剔除unclassified

results1[, 'sum'] <-  rowSums(results1[, mapping$Sample])
results1 <- results1[order(-results1$sum), ]
results1 <- results1[c(which(results1$pvalue <0.05), which(results1$pvalue >=0.05)), ]

if (nrow(results1) == 0){
  print("没有显著差异菌！！！")
  q(0)
}else if(nrow(results1) >30 ){
  results1 <- results1[1:30, ]
}
# results1[, 'sum'] <-  rowSums(results1[, mapping$Sample])
# results1 <- results1[order(-results1$sum), ]
results1[, name] <- factor(results1[, name], levels = results1[, name])
results1 <- results1[, c(mapping$Sample, 'pvalue', "1", name)]
results1[, 'plabel'] <- paste0("p = ", round(results1$pvalue, 4))
results2 <- reshape2::melt(results1, measure=mapping$Sample)
results2 <- merge(results2, mapping, by.x="variable", by.y = "Sample")
colnames(results2) <- c("Sample", "pvalue", "Kingdom", "id", "plabel", "Abundance", "Group")
a <- results2
a <- a[a$Abundance!=0,]
#2.绘图

# newname<-gsub(".txt","",args[1])
# newname<-gsub(".xls","",newname)
# newname<-gsub(".*/","",newname)
# a<-read.table(args[1],sep="\t",header=T)
# gp<-read.table(args[3],sep="\t",header=T)
# colnames(gp) <- c("Sample","Group")
# a['axis']<-paste('p-value',format(a$pvalue,digits=3,scientific=T),sep='=')
gp <- mapping

a$p.value<-as.factor(a$pvalue) 
fa<- unique(gp$Group)
a$Group<-factor(a$Group,levels=fa)
fa2<- unique(a$id)
a$id<- factor(a$id,levels=fa2)

# # 31colors_excel
# mycol <- c("#0B5FD1","#C83406","#008039","#226ED4","#CC4925","#378E4C","#397DD8","#D05F45","#6F9C5F","#508CDB","#D47564","#A6AA72","#679BDF","#D88B84","#DEB885","#7FAAE2","#DB9E9D","#FAC397","#96B9E6","#DDAEB1","#FBCDA7","#ADC8E9","#DFBEC4","#FBD6B8","#C4D7ED","#E1CED8","#FCDFC8","#DCE6F1","#E4DFEC","#FDE9D9","#3466A8")
# singler
mycol <- c("#7FC97F","#BEAED4","#FDC086","#386CB0","#F0027F","#BF5B17",
           "#666666","#1B9E77","#7570B3","#66A61E", "#E6AB02","#A6761D",
           "#A6CEE3","#B2DF8A","#FB9A99","#E31A1C","#FF7F00","#6A3D9A",
           "#8DA0CB","#4DAF4A","#984EA3","#c6c386","#999999","#66C2A5",
           "#FC8D62","#A6D854","#FFD92F","#BEBADA","#FB8072","#80B1D3",
           "#FDB462","#BC80BD","#B3B3B3","#33A02C","#B3DE69","#4038b0",
           "#ee7576","#e94749","#E78AC3","#ff0000","#A65628","#d80172",
           "#F781BF","#D95F02","#E7298A","#1F78B4","#FDBF6F","#CAB2D6",
           "#B15928","#FBB4AE", "#B3CDE3",'#0173b2','#de8f05','#029e73',
           '#d55e00','#cc78bc','#ca9161','#fbafe4','#949494','#ece133',
           '#56b4e9',"#00AFBB","#E7B800","#FC4E07","#FFDB6D","#C4961A",
           "#F4EDCA","#D16103","#C3D7A4","#52854C","#4E84C4","#293352")
mycol2 <- c("#0072B5", "#BC3C29")
# # cold
# mycol -> c("#1660A7","#FF6A00","#219418","#CD0C18","#814BB2","#794339",
#          "#5ddb53","#5066d1","#FF0000","#F0E442","#4836be","#EE6677",
#          "#DC59B6","#CC79A7","#982f29","#AFB400","#11B3C6","#292344",
#          "#E69F00","#AA3377","#009E73","#576e8b","#0072B2","#D55E00",
#          "#00AFBB","#00FFFF","#4477AA","#228833","#CCBB44","#66CCEE",
#           "#56B4E9","#BBBBBB")
my_levels <- levels(a$Group)

# ggplot(data=a, aes(x=id,y=log(Abundance * 100, 10)))+geom_boxplot(aes(fill=Group),outlier.size = 0)+
# #ylim(0,0.01)+
# #facet_wrap( ~ id, scales="free",nrow=1)+
# scale_fill_manual(values=mycol[1:length(my_levels)])+
# labs(title=newname,x=args[2],y="lg(Relative Abundance * 100)")+theme_bw()+
# theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=10,angle = 55,vjust=1,hjust=1))
# ggsave(args[4],width=m,limitsize=FALSE)

####
teda <- log(a$Abundance * 100, 10)
teda[is.infinite(teda)] <- 0
x1 <- max(teda)
x2 <- min(teda)
x <- x1 - x2

y <- length(unique(a$id))

ratio.display <- (7.5/10*y)/2.7
ratio.values <- y/x
ratio_n <- ratio.values / ratio.display

ymin <- min(log(a$Abundance * 100, 10)[log(a$Abundance * 100, 10) != -Inf])
ymax <- max(log(a$Abundance * 100, 10)[log(a$Abundance * 100, 10) != -Inf])
unite_p_y <- (ymax-ymin)/100
ymin <- ymin-unite_p_y*68
ymax <- ymax+unite_p_y*68

p <- ggplot(data = a,mapping = aes(x=id,y=log(Abundance * 100, 10),
                                   colour=Group), fill=NULL)+
  # p <- ggplot(data = a,mapping = aes(x=id,y=Abundance, 
  #                                    colour=Group), fill=NULL)+
  stat_boxplot(mapping=aes(colour=Group),
               size = 1.5,  #定义线大小
               geom = "errorbar",  #添加误差线
               width=0.4,  #设置宽度
               #color="black",  #设置颜色
               position = position_dodge(0.65))+  #调整位置
  geom_boxplot(
    fill="white",
    alpha = 1,  #定义透明度
    size = 1.5,  #定义线的宽度
    width = 0.55,  #定义箱体的宽度
    fatten = 1.8,  #定义中位线的宽度
    outlier.shape = NA,  #不显示异常值
    position = position_dodge(0.65))+  ##调整位置
  geom_jitter(aes(fill=Group),shape = 16,size=1.2, 
              position = position_jitterdodge(jitter.width = 0.5,
                                              dodge.width = 0.75), #width =0.25,
              show.legend = F)+ ##添加散点
  #xlab(args[2])+  ##x轴标题
  xlab("")+
  ylab("lg(Abundance * 100)")+
  # ylim(ymin, ymax)+
  # labs(title = newname)+
  #labs(title = )+
  scale_y_continuous(expand = c(0,0),
                     breaks=seq(round(x2-x/10,digits=2), round(x1+x/10,digits=2), round(x/4,digits=2)),limits=c(x2-x/10, x1+x/10))+ ##设置y轴刻度
  theme_classic(base_line_size = 1)+  #定义绘图主题和基础线条大小
  theme(text=element_text(family="serif"),  #设置字体
        plot.title = element_text(size = 19,  #设置标题字体大小
                                  colour = "black",  #设置标题颜色
                                  face = "bold",  #标题加粗
                                  hjust = 0.5),  #设置标题居中
        
        axis.title.y = element_text(size = 20,  #设置y轴字体大小
                                    color = "black",  #设置y轴字体颜色
                                    #face = "bold",  #设置y轴标题加粗
                                    vjust = 1.9,  #调整y轴标题和轴线的距离
                                    hjust = 0.5,  #设置y轴标题居中
                                    angle = 90),  #设置y轴标题旋转角度
        
        axis.title.x = element_text(size = 19,  #设置y轴字体大小
                                    color = "black",  #设置y轴字体颜色
                                    #face = "bold",  #设置y轴标题加粗
                                    vjust = 1.9,  #调整y轴标题和轴线的距离
                                    hjust = 0.5,  #设置y轴标题居中
                                    angle = 0),  #设置y轴标题旋转角度
        
        legend.title = element_text(color="black",  #修改图例的标题颜色
                                    size=16,  #修改图例的标题大小
                                    face="bold"),  #图例的标题加粗
        
        legend.text = element_text(color="black",  #设置图例标签文字
                                   #face="bold",
                                   size = 20),  #设置图例标签字体大小
        legend.position = c(1.05,0.7),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = 15,  #设置X轴标签字体大小
        #                            color = "black",  #设置X轴标签颜色
        #                            vjust = 1,  #设置X轴标签与轴线距离
        #                            hjust = 1,  #设置X轴标签居中
        #                            angle = 45,  #设置X轴标签旋转角度
        #                            face="bold"),
        
        axis.text.y = element_text(size = 18,  #设置y轴标签字体大小
                                   color = "black",  #设置y轴标签颜色
                                   vjust = 0.5,  #设置y轴标签与轴线距离
                                   hjust = 0.5,  #设置y轴标签居中
                                   angle = 0),  #设置y轴标签旋转角度
        
        axis.line.x = element_line(size=1),  #设置X轴轴线大小
        axis.line.y = element_line(size=1),  #设置X轴轴线大小
        axis.ticks.x= element_line(size=1),  #设置X轴刻度线大小
        axis.ticks.y= element_line(size=1))+  ##设置y轴刻度线大小
  theme(plot.margin = margin(t = 0.5,  # 顶部边缘距离
                             r = 5,  # 右边边缘距离
                             b = 8,  # 底部边缘距离
                             l = 3,  # 左边边缘距离
                             unit = "cm") #设置单位为cm
  )+
  scale_fill_manual(values=mycol2[1:length(my_levels)]) + 
  scale_colour_manual(values=mycol2[1:length(my_levels)]) + 
  theme(legend.title=element_blank())+
  # coord_fixed(ratio = ratio_n)+ #宽高比
  coord_cartesian( expand = T, clip = "off")+
  scale_y_continuous(limits = c(ymin, ymax))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))  ###设置画框

#文字添加颜色
col = data.frame("Kingdom" = c("Archaea", "Bacteria", "Fungi", "Viruses"),
                 color = mycol[1:4])
col <- merge(unique(a[, c("Kingdom", "id")]), col, by.x = "Kingdom", by.y = "Kingdom")
row.names(col) <- col$id

library(gridExtra)
library(grid)
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}
for(i in levels(a$id)){
  p <- p+annotation_custom2(
    grob = textGrob(i, 
                    gp=gpar(col=col[i, 'color'], lwd=1, cex = 2,
                            lty="dotted", fontfamily="serif"#,fontface="bold"
                            ),
                    rot=45,
                    hjust=1, vjust=1), 
    data = a,
    xmin = i, xmax = i, 
    ymin = ymin-unite_p_y*15, ymax = ymin-unite_p_y*15)
}

#添加圖例

col[, 'y'] <- 1:nrow(col)
p_tmp <- ggplot(col, aes(x=id, y=y, color=Kingdom)) +
  geom_point(shape=15)+ scale_color_manual(values=unique(col$color), name="")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(text=element_text(family="serif"),  #设置字体
        legend.title = element_text(color="black",  #修改图例的标题颜色
                                    size=16,  #修改图例的标题大小
                                    face="bold"),  #图例的标题加粗
        legend.box.background = element_blank(),
        # legend.key.size = unit(5,"cm"),
        legend.key = element_blank(),
        legend.text = element_text(color="black",  #设置图例标签文字
                                   #face="bold",
                                   size = 20)  #设置图例标签字体大小
  )
legend_t <- cowplot::get_legend(p_tmp)
p <- p + ggimage::geom_subview(x=nrow(col)*1.074, y=ymin+unite_p_y*70, subview=legend_t)

#添加P值
library(dplyr)

up_degree <- 36
floor_degree <- 36
a[, 'Abundance_log'] <- log(a$Abundance * 100, 10)
p_text <- a %>% group_by(id) %>% 
  summarise(p_text_y = if((ymin+ymax)/2 > median(Abundance_log)) {
    unite_p_y*up_degree+max(Abundance_log)  #p值在箱子上
  }else{
    min(Abundance_log)-unite_p_y*floor_degree  #p值在箱子下
  },
  median = median(Abundance_log),
  max = max(Abundance_log),
  min = min(Abundance_log))

results1[, 'plabel'] <- paste0("p = ", round(results1$pvalue, 4))
p_text <- merge(p_text, results1[,c(name, 'plabel')], by.x='id', by.y=name)
p_text$id <- factor(p_text$id, levels = levels(a$id))
p <- p+geom_text(data =p_text, aes(x=id, y=p_text_y,
                              label = plabel ),alpha=1,
                 fontface="italic", color="grey",angle=90,size=4
)


####保存图片
gp <- mapping
nu = length(unique(a$id))
group = length(unique(a$Group))

margin_l = max((nchar(levels(a$id))-33)*0.2+0.5, 0.5)
if (nu>5){
  m=nu*group*0.5
}else if (nu>1 & group>4){
  m=nu*group*0.5
}else{
  m=7
}

m = group*2 + 12

pdf(paste0(args[3], "/L",name ,".diff.top30.boxplot.pdf"), width = 20.4, height = 7.4)
print(p)
dev.off()


# write.table(heatmap,args[4],sep="\t",quote=F,row.names=F)
# 
# outdir <- dirname(gsub(".diff.heatmap.xls", "", args[4]))
# outname <- basename(gsub(".diff.heatmap.xls", "", args[4]))
# 
# if(args[5] == "w" || args[5] == "t") {
#   s1 <- rownames(mapping)[
#     mapping[, 1, drop = T] == unique(
#       mapping[2:nrow(mapping), "Group", drop = T])[1]
#   ]
#   s2 <- rownames(mapping)[
#     mapping[, 1, drop = T] == unique(
#       mapping[2:nrow(mapping), "Group", drop = T])[2]
#   ]
#   m1 <- rowSums(heatmap[, s1, drop = F]) / length(s1)
#   m2 <- rowSums(heatmap[, s2, drop = F]) / length(s2)
#   lfc <- log2(m1 / m2)
#   omics <- merge(heatmap, results[,c('Index', 'pvalue'),drop=F], by.x='Index',by.y='Index')
#   omics$lfc <- lfc
#   write.table(
#     omics,
#     paste(
#       c(
#         file.path(outdir, gsub("\\..*", "", outname)),
#         paste(unique(mapping[2:nrow(mapping), "Group", drop = T]), collapse = "-"),
#         "txt"
#       ),
#       collapse = "."
#     ),
#     sep = "\t", quote = F, row.names = F
#   )
# }