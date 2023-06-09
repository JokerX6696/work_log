rm(list=ls())
setwd('D:\\desk\\DZOE2023041442-王石-数据高级分析')
library(ggplot2)
library(reshape2)
df <- read.table('stat_group_2.txt',sep = '\t',header = TRUE)
df <- df[-1,]
df$sample <- c('Cyprinus carpio','Megalobrama amblycephala')
data <- melt(df,id.vars = 'sample')
data$variable <- gsub('\\.', '-', data$variable)
colnames(data)[1] <- 'Species'
COL <- c('#FC8C64','#8EA0CB')

ggplot(data = data,mapping = aes(x=factor(variable,levels = rev(unique(data$variable))), y = value, fill = Species)) + 
  geom_bar(stat = 'identity',colour="#414141",width=0.6,position = "fill") +
  labs(x="",y="Percent") +
  scale_fill_manual(breaks = data$Species,values = COL) +
  scale_y_continuous(breaks=c(0.00,0.25,0.50,0.75,1.00), labels = c('0%','25%','50%','75%','100%')) +
  #geom_text(data = data,mapping = aes(x=factor(variable,levels = rev(unique(data$variable))), y = value),label = 2) +
  coord_flip() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5) )
ggsave(filename = 'All_sample_2_group.png',device = 'png',width = 8, height = 6)
ggsave(filename = 'All_sample_2_group.pdf',device = 'pdf',width = 8, height = 6)
 