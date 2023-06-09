args<-commandArgs(T)

# library(hash)
library(ggplot2)
library("dplyr")
library("ggsignif")
library(plotrix)

# setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
# args <- c("QC.xls", "mapping.txt", "clean_bases")

data <- read.table(args[1],header=T,sep="\t",row.names=1,check.names=FALSE)
group = read.table(args[2], sep="\t", header=T, check.name = FALSE, quote = "", comment.char ="")
# group <- group[, c(1,4)]
colnames(group) = c("Samples","Group")
#做哪几个组
# l = c('A', 'B', 'N')
# group <- group[group$Group %in% l, ]
# data <- data[,colnames(data) %in% group$Samples]
# 
# #组名
# lab1 <- c('A', 'B', 'N')
# group1 <- hash(keys= lab1, values=c("T", "NA", "N"))
# #组名
# lab <- c()
# for (j in unique(group$Group)) {
#   lab <- c(lab, group1[[j]])
# }
# #颜色
# board <- hash(keys=lab1, values=c("#038656", "#8C6B03", "#005CAE"))
# #颜色
# board1 <- c()
# for (j in unique(group$Group)) {
#   board1 <- c(board1, board[[j]])
# }
# #输出路径
# output <- paste0(unique(group$Group), collapse= "_")
# 
#图片尺寸
group_num = length(levels(group$Group))
if (group_num < 3){
  width = 6.8 + 1.3*group_num
}else{
  width = 6.8 + 1.0*group_num
}
height = 7.4 + 0.2*group_num


library(stringr)
Clean_Base1 <- c()
for (i in rownames(data)){
  Clean_Base1 <- c(Clean_Base1, str_match(data[i, "Clean_Base"], "\\((.*)G\\)")[,2])
}
Clean_Base1 <- as.numeric(Clean_Base1)
df = data.frame("Clean_Base"=Clean_Base1, "Samples" = row.names(data))
plotdata <- merge(df, group, by.x = 'Samples', by.y = 'Samples')
colnames(plotdata) <- c("sample", "abundance", "Group")
plotdata$Group<-factor(plotdata$Group,levels=unique(group$Group))
# plotdata <- plotdata[plotdata$abundance != 0,]

tongji <- plotdata %>% group_by(Group) %>% summarise(value = sum(abundance))
plotdata <- plotdata[!(plotdata$Group %in% tongji$Group[tongji$value %in% 0]), ]



res <- pairwise.wilcox.test(plotdata$abundance, plotdata$Group, p.adjust.method='none')
# res <- pairwise.t.test(plotdata$abundance, plotdata$Group, p.adjust.method='none')
output <- args[3]
if(!file.exists(output)){dir.create(output)}
write.table(res$p.value, paste0(output,"/Clean_Bases_bar.xls"),quote=F,sep="\t",col.names = NA)


pValue <- c(res$p.value)[!is.na(c(res$p.value))]
my_annotations <- as.character(symnum(pValue, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***",'**','*','')))
my_annotations1 <- my_annotations[my_annotations != ""]

wtq <- levels(plotdata$Group)
lis <- combn(wtq,2)
my_comparisons <- tapply(lis, rep(1:ncol(lis), each=nrow(lis)), function(i) i)
my_comparisons1 <- my_comparisons[my_annotations != ""]

#bar
# library(tidyverse)
# detach('dplyr')
bar <-plotdata %>%
  group_by(Group) %>%
  summarise(
    median = median(abundance),
    sd = sd(abundance),
    std_e = std.error(abundance),
    ymax = max(abundance),
    ymin = min(abundance),
    yq25 = quantile(abundance,0.25),
    yq75 = quantile(abundance,0.75),
    abundance  = mean(abundance))
bar$Group <- factor(bar$Group, levels = bar$Group)
bar$id <- 1:length(bar$Group)
bar$xmin <- bar$id-0.3
bar$xmax <- bar$id+0.3
#颜色
# board1 <- c()
# for (j in unique(bar$Group)) {
#   board1 <- c(board1, board[[j]])
# }
board1 <- c("#35905F", "#E3923C")
#组名
lab <- c()
# for (j in unique(bar$Group)) {
#   lab <- c(lab, group1[[j]])
# }
p <- ggplot(data=bar, aes(x=Group,y=abundance,fill=Group))+ 
  geom_bar(colour="black",size=1.2,alpha = 1, width = 0.7, stat = "identity", show.legend=F)+
  geom_errorbar(aes(ymin=abundance-std_e, ymax=abundance+std_e), 
                width=0.4, alpha=1, color="black",
                size = 1.2, show.legend=F) +
  # geom_jitter(alpha = 1, size = 2, width = 0.25, height = 0, stat = "identity")+
  # geom_segment(data=bar, aes(x=xmin, xend=xmax,  y=median, yend=median, group=Group, alpha=0.5), 
  #              size=1.5, lineend="round", linejoin="round", show.legend=F)+
  # geom_segment(data=bar, aes(x=xmin, xend=xmax,  y=yq25, yend=yq25, group=Group, alpha=0.5), 
  #              size=1.5, lineend="round", linejoin="round", show.legend=F)+
  # geom_segment(data=bar, aes(x=xmin, xend=xmax,  y=yq75, yend=yq75, group=Group, alpha=0.5), 
  #              size=1.5, lineend="round", linejoin="round", show.legend=F)+
  labs(title="", x="", y="Bases (G) per sample")+
  # scale_x_discrete(labels = lab)+
  scale_fill_manual(labels = levels(bar$Group), values = board1)+
  theme_bw()+
  coord_cartesian(clip = 'off')+
  theme(text = element_text(family = "serif"),
    plot.title = element_text(size = 22, face = "bold.italic",hjust=0.5), 
        axis.title.x = element_text(size = 0, face = "bold"), 
        axis.title.y = element_text(size = 40), 
        axis.text.x = element_text(size=40,angle = 0, colour = "black", vjust=0.5, hjust=0.5),
        axis.text.y = element_text(size=36,angle = 0, colour = "black", vjust=0.5, hjust=0.5),
        axis.ticks =element_line(size = 1.2) ,
        axis.ticks.x.top = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line = element_line(colour ="black", size = 1.2), 
        axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.x.top  = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black",size = 1.2),
        axis.line.y.right = element_line(colour = "black", size = 1.2),
        legend.title=element_blank(), 
        legend.text=element_text(size=18), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid = element_blank(),
        # panel.background = element_blank())
        panel.background = element_rect(colour = "black", size = 1.5))

if (length(my_annotations1) > 0) p = p + geom_signif(annotations = my_annotations1, comparisons = my_comparisons1, step_increase =0.1, vjust=0.7, colour='gray20', margin_top = 0.3,tip_length=0.03, size=1, textsize= 10)


pdf(file=paste0(output,"/Clean_Bases_bar.pdf"), width=width, height=height)
print(p)
dev.off()

