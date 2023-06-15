rm(list=ls())
setwd('D:\\desk\\XMSH_202305_4038')
args<- c("otu_table_even_L2.txt",
         "mapping.txt",
         "diff_boxplot",
         "t")
library(dplyr)

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

HCC <- HCC[,-c(1:6)]
HCC_df <- data.frame(Specie=HCC$Specie,HCC=apply(HCC[,2:length(HCC)], 1,mean))
cancer <- cancer[,-c(1:6)]
cancer_df <- data.frame(Specie=cancer$Specie,Cancer=apply(cancer[,2:length(cancer)], 1,mean))


df <- merge(x=HCC_df,y=cancer_df,by='Specie',all = F)


rownames(df) <- df$Specie

df <- df[,-1]
df <- df[rev(order(apply(df, 1, sum))),]
df <- df[1:30,]

library(pheatmap)
annotation_col = data.frame(Group = factor(c("HCC", "Cancer")))
rownames(annotation_col) = c("HCC", "Cancer")
pdf(file = '2bM_heatmap.pdf',height = 9,width = 6)
pheatmap(df,
         #scale = "row",
         angle_col = 0,
         #cellwidth = 15,
         #fontsize_row=2,
         show_colnames=T,
         annotation_col = annotation_col,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)

dev.off()









# 
# 
# #########################
# taxon <- c()
# for(i in 1:dim(HCC)[1]){
#   t <- paste(
#     paste('k',HCC[i,1],sep = '__'),
#     paste('p',HCC[i,2],sep = '__'),
#     paste('c',HCC[i,3],sep = '__'),
#     paste('o',HCC[i,4],sep = '__'),
#     paste('f',HCC[i,5],sep = '__'),
#     paste('g',HCC[i,6],sep = '__'),
#     paste('s',HCC[i,7],sep = '__'),
#     sep = ';'
#   )
#   taxon <- c(taxon,t)
# }
# HCC <- data.frame(Taxon=taxon,HCC[8:length(HCC)])
# LEVEL = 'Taxon'
# #HCC <- fun(HCC,LEVEL)
# rownames(HCC) <- HCC$Taxon 
# 
# taxon <- c()
# for(i in 1:dim(cancer)[1]){
#   t <- paste(
#     paste('k',cancer[i,1],sep = '__'),
#     paste('p',cancer[i,2],sep = '__'),
#     paste('c',cancer[i,3],sep = '__'),
#     paste('o',cancer[i,4],sep = '__'),
#     paste('f',cancer[i,5],sep = '__'),
#     paste('g',cancer[i,6],sep = '__'),
#     paste('s',cancer[i,7],sep = '__'),
#     sep = ';'
#   )
#   taxon <- c(taxon,t)
# }
# cancer <- data.frame(Taxon=taxon,cancer[8:length(HCC)])
# #cancer <- fun(cancer,LEVEL)
# rownames(cancer) <- cancer$Taxon 
# 
# 
# filedata <- merge(x=cancer,y=HCC,by = 'Taxon')
# rownames(filedata) <- filedata$Taxon
# filedata <- filedata[,-1]
# 
# mapping <- data.frame(Sample=names(filedata),Group=sub("_.+","",names(filedata)))
# rownames(mapping) <- names(filedata)
# # ####################
# # LEVEL = q
# # HCC <- fun(HCC,LEVEL)
# # cancer <- fun(cancer,LEVEL)
# # HCC <- HCC[apply(HCC[,-1],1,sum) != 0,]
# # cancer <- cancer[apply(cancer[,-1],1,sum) != 0,]
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# ########################################################################
