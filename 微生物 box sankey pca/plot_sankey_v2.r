rm(list=ls())
setwd('D:\\desk\\XMSH_202305_4038')
args<- c("otu_table_even_L2.txt",
         "mapping.txt",
         "diff_boxplot",
         "t")
library(dplyr)
l=3
k=1.22
q='Phylum'
fun <- function(df,LEVEL){
  df_sum <- df %>%
    group_by_at(vars(LEVEL)) %>%
    summarize(across(names(df)[8:dim(df)[2]], sum, .names = "{.col}"))
  return(df_sum)
}
#for(q in c('Phylum','Genus','Specie')){}
# HCC 处理
HCC <- read.table('宏基因组Abundance_Stat.txt',quote = "", sep = '\t',header = T)
HCC <- HCC[,c(1:7,grep('HCC',names(HCC)))]
for(num in seq(8,dim(HCC)[2]) ){
  names(HCC)[num] = paste('HCC','_',num-7,sep = '')
}


# cancer 处理
cancer <- read.table('Abundance_Stat.filter.anno.xls',quote = "", sep = '\t',header = T)
cancer <- cancer[,c(1:7,grep('Tumor|HCC',names(cancer)))]
for(num in seq(8,dim(cancer)[2]) ){
  names(cancer)[num] = paste('cancer','_',num-7,sep = '')
}





#########################
taxon <- c()
for(i in 1:dim(HCC)[1]){
  t <- paste(
    paste('k',HCC[i,1],sep = '__'),
    paste('p',HCC[i,2],sep = '__'),
    paste('c',HCC[i,3],sep = '__'),
    paste('o',HCC[i,4],sep = '__'),
    paste('f',HCC[i,5],sep = '__'),
    paste('g',HCC[i,6],sep = '__'),
    paste('s',HCC[i,7],sep = '__'),
    sep = ';'
  )
  taxon <- c(taxon,t)
}
HCC <- data.frame(Taxon=taxon,HCC[8:length(HCC)])
LEVEL = 'Taxon'
#HCC <- fun(HCC,LEVEL)
rownames(HCC) <- HCC$Taxon 

taxon <- c()
for(i in 1:dim(cancer)[1]){
  t <- paste(
    paste('k',cancer[i,1],sep = '__'),
    paste('p',cancer[i,2],sep = '__'),
    paste('c',cancer[i,3],sep = '__'),
    paste('o',cancer[i,4],sep = '__'),
    paste('f',cancer[i,5],sep = '__'),
    paste('g',cancer[i,6],sep = '__'),
    paste('s',cancer[i,7],sep = '__'),
    sep = ';'
  )
  taxon <- c(taxon,t)
}
cancer <- data.frame(Taxon=taxon,cancer[8:length(HCC)])
#cancer <- fun(cancer,LEVEL)
rownames(cancer) <- cancer$Taxon 


filedata <- merge(x=cancer,y=HCC,by = 'Taxon')
rownames(filedata) <- filedata$Taxon
filedata <- filedata[,-1]

mapping <- data.frame(Sample=names(filedata),Group=sub("_.+","",names(filedata)))
rownames(mapping) <- names(filedata)
# ####################
# LEVEL = q
# HCC <- fun(HCC,LEVEL)
# cancer <- fun(cancer,LEVEL)
# HCC <- HCC[apply(HCC[,-1],1,sum) != 0,]
# cancer <- cancer[apply(cancer[,-1],1,sum) != 0,]
# 








########################################################################
row.names(mapping) <- mapping$Sample
num <- length(rownames(filedata))


#1.计算P值

res <- NULL
for(i in 1:num){
  if(sd((filedata[i,]))!=0){
    x <- as.data.frame(t(filedata[i,]))
    colnames(x) <- c("datainfo")
    merged_data<-merge(x,mapping[,'Group',drop=F],by.x='row.names',by.y='row.names')
    if('t'=="aov"){
      method <- "ANOVA"
      merged_data.aov <- aov(datainfo ~ Group,data=merged_data)
      pvalue <- summary(merged_data.aov)[[1]][5][1,]
    }else if('t'=="kw"){
      method <- "Kruskal Wallis"
      p <- kruskal.test(datainfo ~ Group,data=merged_data)
      pvalue <- p$p.value
    }else if('t'=="t"){
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
    }else if('t'=="w"){
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
#差异 top30
results1 <- res[res$pvalue < 0.05,]
results1 <- results1[!grepl("s__uniclassified", row.names(results1)), ]
results1 <- results1[!grepl("g__uniclassified", row.names(results1)), ]
results1 <- results1[!grepl("p__uniclassified", row.names(results1)), ]
results1[, 'sum'] <-  rowSums(results1[, mapping$Sample])
results1 <- results1[order(-results1$sum), ]
#top30 <- results1[1:30, ]
top30 <- results1[1:26, ]
top30[, 'mean'] <-  rowMeans(top30[, mapping$Sample])
top30 <- cbind(top30, do.call(rbind, strsplit(row.names(top30), ";")))
for(i in 2:7){
  top30[, as.character(i)] <- do.call(rbind, strsplit(top30[, as.character(i)], "__"))[,2]
  
}
top30 <- top30[, c('2', '6', '7', 'mean')]
colnames(top30) <- c('Phylum', 'Genus', 'Species', 'mean')

#term分组
group_res = data.frame()
for(i in c('Phylum', 'Genus', 'Species')){
  if(nrow(group_res) == 0){
    group_res = data.frame("env"=unique(top30[,i]), 'group'=i)
  }else{
    group_res = rbind(group_res, data.frame("env"=unique(top30[,i]), 'group'=i))
  }
  
}
group_res[, 'description'] <- 1:nrow(group_res)-1
#

df <- top30 %>% group_by(Genus) %>%
  summarise(value = sum(mean),
            name = unique(Phylum)
  )
colnames(df) <- c('source','weight','target' )
df1 <- top30[, c('Species', 'Genus', 'mean')]
colnames(df1) <- c('source','target','weight')
df <- rbind(df[,c('source','target','weight')], df1)
df <- merge(df, group_res, by.x='source', by.y='env')
colnames(df) <- c('source','target','weight', 'source_g', 'IDsource')
df <- merge(df, group_res, by.x='target', by.y='env')
colnames(df) <- c('target','source','weight', 'source_g', 'IDsource', 'target_g', 'IDtarget')


#
# #“env_table.txt” 记录了测量的环境变量数据
# env <- read.delim(args[1], row.name = 1, check.names = FALSE)
# env <- t(env)
#
# #量纲不同，作个标准化，尽管标准化前后对相关系数的计算无影响
# env <- scale(env)
#
# #环境变量间的相关性，以 spearman 秩相关系数为例
# env_corr <- rcorr(t(env), type = 'spearman')
#
# r <- env_corr$r  #相关系数 r 值
# p <- env_corr$P  #相关系数的显著性 p 值
# p <- p.adjust(p, method = 'BH')  #可选 p 值校正，这里使用 BH 法校正 p 值
# p1 <- matrix(p, nrow=nrow(r), ncol=ncol(r))
# rownames(p1) <- rownames(r)
# colnames(p1) <- colnames(r)
# #输出相关系数矩阵
# write.table(r, 'corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
# write.table(p1, 'pvalue.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
#
# # #仅保留 r>=0.7 且 p.adj<0.01 的相关系数，且去除对角线的自相关
# # r[abs(r)<0.7 | p>=0.01] <- 0
# # diag(r) <- 0
# #仅保留p.adj<0.05 的相关系数，且去除对角线的自相关
# # r[p >= 0.05] <- 0
# #
# # #输出相关系数矩阵
# # write.table(r, 'env_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
#
# # #使用 corrplot 包绘制相关图可视化相关矩阵
# # library(corrplot)
# # col1 <- colorRampPalette(c('blue4', 'blue', 'white', 'orange', 'red3'))
# #
# # corrplot(r, method = 'number', col = col1(21),
# #          tl.cex = 0.8, number.cex = 0.8, cl.length = 5, diag = FALSE)
# # corrplot(r, method = 'pie',col = col1(21),
# #          tl.pos = 'n', cl.pos = 'n', add = TRUE, type = 'upper', diag = FALSE)
# #
#
# #################################
# #        2.构建桑基图数据
# #################################
#
# #读取刚才计算的相关矩阵
#
# r <- read.delim('corr.matrix.txt', check.names = FALSE)
# p <- read.delim('pvalue.matrix.txt', check.names = FALSE)
#
# #由于这是一个对称矩阵，仅保留下三角数值即可，把上三角的冗余值去除
# r[upper.tri(r)] <- 0
# # r[p >= 0.05] <- 0
# #将相关系数矩阵转为“两两关系列表”的样式
# env_corr_list <- reshape2::melt(r)
# names(env_corr_list) <- c('source', 'target', 'spearman')
# env_corr_list$source <- as.character(env_corr_list$source)
# env_corr_list$target <- as.character(env_corr_list$target)
# env_corr_list <- env_corr_list[env_corr_list$spearman != 0,]
#
# #添加P值
# row.names(p) <- p[,1]
# for(i in 1:nrow(env_corr_list)){
#   env_corr_list[i, "pvalue"] <- p[env_corr_list[i,"source"], env_corr_list[i,"target"]]
# }
#
# #剔除同组的关联
# env_group <- read.delim('env_group.txt', stringsAsFactors = FALSE, check.names = FALSE)
# env_corr_list <- merge(env_corr_list, env_group[,1:2], by.x = "source", by.y ="env")
# env_corr_list <- rename(env_corr_list, source_g=group)
# env_corr_list <- merge(env_corr_list, env_group[,1:2], by.x = "target", by.y ="env")
# env_corr_list <- rename(env_corr_list, target_g=group)
# env_corr_list <- env_corr_list[env_corr_list$source_g != env_corr_list$target_g,]
#
# #剔除KO与Behavioral_Indicators的关联
# env_corr_list <- env_corr_list[
#   !(env_corr_list$source_g =="KO" & env_corr_list$target_g=="Behavioral_Indicators") &
#     !(env_corr_list$source_g =="Behavioral_Indicators" & env_corr_list$target_g=="KO"),]
#
# #KO放在左边，代谢放在中间，Behavioral_Indicators放在右边
# tmp1 <- env_corr_list[(env_corr_list$source_g =="Behavioral_Indicators" & env_corr_list$target_g=="Metabolites"),]
# tmp2 <- env_corr_list[!(env_corr_list$source_g =="Behavioral_Indicators" & env_corr_list$target_g=="Metabolites"),]
# tmp1 <- rename(tmp1, target=source, source=target, target_g=source_g, source_g=target_g)
# env_corr_list <- rbind(tmp1, tmp2)
#
#
#
# #取相关系数的绝对值赋值新列作为相关的强度，并取相关系数的符号赋值新列作为相关的方向
# #以及去除 0 值
# env_corr_list$weight <- abs(env_corr_list$spearman)
# env_corr_list[which(env_corr_list$spearman >= 0.8 & env_corr_list$pvalue < 0.05),'direction'] <- '1'
# env_corr_list[which(env_corr_list$spearman <= -0.8 & env_corr_list$pvalue < 0.05),'direction'] <- '-1'
# env_corr_list[env_corr_list$direction %in% NA,'direction'] <- '0'
# env_corr_list$label <- as.character(round(env_corr_list$spearman, 3))
#
# #读取环境变量的属性列表，分配 id 指代环境变量，并用在后续根据分类赋值颜色等
# env_group <- subset(env_group, env %in% unique(c(env_corr_list$source, env_corr_list$target)))
# env_corr_list$IDsource <- match(env_corr_list$source, env_group$env) - 1
# env_corr_list$IDtarget <- match(env_corr_list$target, env_group$env) - 1
#
# head(env_corr_list)
#
# #输出相关性列表
# write.table(env_corr_list, 'env_corr.sankeynetwork.txt', row.names = FALSE, sep = '\t', quote = FALSE)
# write.table(env_group, 'env_group.sankeynetwork.txt', row.names = FALSE, sep = '\t', quote = FALSE)
#
# ###########################################
# #           3.绘图
# ###########################################
# ##3.1 networkD3 包的桑基图

# #读入绘图文件
# #读入输入文件
# env_corr_list <- read.delim('env_corr.sankeynetwork.txt', check.names = FALSE)
# env_group  <- read.delim('env_group.sankeynetwork.txt', check.names = FALSE)

env_corr_list <- df
env_group  <- group_res
#标签和方向必须为字符串
# env_corr_list$direction <- as.character(env_corr_list$direction)
# env_corr_list$label <- as.character(env_corr_list$label)

#排序
# env_corr_list <- env_corr_list[order(env_corr_list$source_g,
#                                      env_corr_list$IDsource,
#                                      env_corr_list$target_g,
#                                      env_corr_list$IDtarget), ]
# env_corr_list <- env_corr_list[order(env_corr_list$source_g), ]
# env_corr_list$source_g <- factor(env_corr_list$source_g, levels = c('Phylum', 'Genus', 'Species'))
data <- env_corr_list
# data <- data.frame()
# for(i in unique(env_group$group)){
#   data <- rbind(data,env_corr_list[env_corr_list$source_g == i,])
# }

#放大线的差异
#v为向量，enlarge为放大到那个范围，为向量；x为输出
FEnlarge <- function(v, enlarge=c(1,15)){
  tmp <-c()
  for(i in v){
    tmp = c(tmp,(enlarge[2]-enlarge[1])/(max(v)-min(v))*i+enlarge[2]-max(v)*(enlarge[2]-enlarge[1])/(max(v)-min(v)))
  }
  return(tmp)
}

# data$weight <- FEnlarge(data$weight)

#定义节点和连线的颜色
#"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#FF00004D", "#0000FF4D"
# c(unique(env_group$group), c("1", "-1"))
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
# color <- 'd3.scaleOrdinal().domain(["KO", "Metabolites", "Behavioral_Indicators", "1", "-1", "0"])
# .range(["#A67D07", "#0818B5", "#0F6A08", "#FABDA4", "#85E8FC", "#B4B2B2"])'

#指定预定义的变量 id，按相关性方向赋值连线颜色，相关性强度赋值连线尺寸
#节点颜色按预定义的变量的分类着色
#节点、字体、边界尺寸等，详情 ?sankeyNetwork
# env_group$group <- factor(env_group$group, levels = c('Phylum', 'Genus', 'Species'))
colnames(data) <-  c('source','target','weight', 'target_g', 'IDtarget', 'source_g', 'IDsource')

if(!dir.exists('sankey_abundance')){
  dir.create('sankey_abundance')
}

write.table(data, paste0('sankey_abundance', '/sankeynetwork.link.xls'), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(env_group, paste0('sankey_abundance', '/sankeynetwork.node.xls'), row.names = FALSE, sep = '\t', quote = FALSE)

p <- sankeyNetwork(Links = data, Nodes = env_group,
                   Source = 'IDsource', Target = 'IDtarget',
                   Value = 'weight', #LinkGroup = 'source_g',
                   NodeID = 'env', NodeGroup = 'group', sinksRight=TRUE,
                   nodePadding = 5, nodeWidth = 25, fontSize = 15,fontFamily ="serif",
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"),
                   margin=2,
                   height =1000, width = 1000)


saveNetwork(p, paste0('sankey_abundance', "/sankey_abundance.html"))

psycModel::html_to_pdf(file_path = paste0('sankey_abundance', "/sankey_abundance.html"),
                       render_exist = T)



