rm(list=ls())
setwd('D:\\desk\\XMSH_202305_4038')
library(dplyr)

fun <- function(df,LEVEL){
  #print(df[,8:dim(df)[2]])
  df_sum <- df %>%
    group_by_at(vars(LEVEL)) %>%
    summarize(across(names(df)[8:dim(df)[2]], sum, .names = "{.col}"))
  return(df_sum)
}
for(q in c('Phylum','Genus','Specie')){
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
  ####################
  LEVEL = q
  HCC <- fun(HCC,LEVEL)
  cancer <- fun(cancer,LEVEL)
  HCC <- HCC[apply(HCC[,-1],1,sum) != 0,]
  cancer <- cancer[apply(cancer[,-1],1,sum) != 0,]
  
  set1 <- cancer[[LEVEL]]
  set2 <- HCC[[LEVEL]]
  
  outname <- paste0(LEVEL,'_venn.pdf')
  library(VennDiagram)
  p <- venn.diagram(
    x = list(A=set1, B=set2),
    category.names = c("Cancer","HCC" ),
    filename = NULL,
    output=FALSE,
    #col = c("#0072B5", "#BC3C29"),
    fill = c("#0072B5", "#BC3C29"),
    #width.prop = 1,
    alpha = 1,
    scaled = FALSE,# 圆圈大小设为固定
    cat.pos = c(12,0),
    cat.cex = 2
  )
  pdf(file = outname,width = 9,height = 9)
  grid.draw(p)
  dev.off()
  
  
}







