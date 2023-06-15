rm(list=ls())
setwd('D:\\desk\\XMSH_202305_4038')
library(grid)
library(ggplot2)
meta <- read.table(file = '宏基因组Abundance_Stat.txt',encoding="UTF-8", sep = "\t", header = TRUE)
toBM <- read.table(file = 'Abundance_Stat.filter.anno.xls', sep = "\t", header = TRUE, encoding = 'UTF-8')
mycol2 <- c("#0072B5", "#BC3C29")
new_df <- merge(x=meta,y = toBM,by = 'Specie',all = F)

no_metch <- toBM[!toBM$Specie %in% new_df$Specie,7]

new_df <- new_df[,-c(64:69)]
# 属性矩阵
lab_df <- data.frame(new_df[2:7],new_df[1])
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
group <- c(rep('cancer',115),rep('paracancer',30))
mapping <- data.frame(Group=group)
rownames(mapping) <- samples
abundence <- pca_df

############################----------------------------------------------------------------
library(pheatmap)
cancer_2 <- apply(cancer,1,mean)
paracancer_2 <- apply(paracancer,1,mean)
heatmap_df <- data.frame(cancer=cancer_2,paracancer=paracancer_2)
heatmap_df$mean <- apply(heatmap_df,1,mean)
rownames(heatmap_df) <- new_df$Specie
heatmap_df <- heatmap_df[rev(order(heatmap_df$mean)),]
heatmap_df <- heatmap_df[1:30,-3]
annotation_col = data.frame(Group = factor(c("cancer", "paracancer")))
rownames(annotation_col) = c("cancer", "paracancer")
pdf(file = 'meta_2bM_heatmap.pdf',height = 9,width = 6)
pheatmap(heatmap_df,
         scale = "row",
         angle_col = 0,
         #cellwidth = 15,
         #fontsize_row=2,
         show_colnames=T,
         annotation_col = annotation_col,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
         )

dev.off()



