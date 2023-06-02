rm(list = ls())
args <- commandArgs(T)

# setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
# args <-c("L7.xls", "mapping.txt", "heatmap_diff_average", "t")
# args <-c("KEGG.level3.abundance.xls", "mapping.txt", "heatmap_KEGG")

library(reshape2)
library(dplyr)
library(ComplexHeatmap)



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
#差异 top30
results1 <- res[res$pvalue < 0.05,]
results1 <- results1[!grepl("s__uniclassified", row.names(results1)), ]
# results1 <- results1[!grepl("g__uniclassified", row.names(results1)), ]
# results1 <- results1[!grepl("p__uniclassified", row.names(results1)), ]
results1[, 'sum'] <-  rowSums(results1[, mapping$Sample])
results1 <- results1[order(-results1$sum), ]
top30 <- results1[1:30, ]
# top30[, 'mean'] <-  rowMeans(top30[, mapping$Sample])
top30 <- cbind(top30, do.call(rbind, strsplit(row.names(top30), ";")))
for(i in 2:7){
  top30[, as.character(i)] <- do.call(rbind, strsplit(top30[, as.character(i)], "__"))[,2]  
  
}
for(i in unique(mapping$Group)){
  top30[, paste0(i)] <-  rowMeans(top30[, mapping[mapping$Group==i, 'Sample']])
}

top30 <- top30[, c('7', unique(mapping$Group))]
colnames(top30) <- c('Species', unique(mapping$Group))



df <- top30
row.names(df) <- df$Species
df <- df[,-1]

if(!dir.exists(args[3])){
  dir.create(args[3])
}

write.table(df, paste0(args[3], '/heatmap_diff_top30_average.xls'), row.names = FALSE, sep = '\t', quote = FALSE)

# 
# 
# 
# 
# mapping<-read.table(args[2],header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=FALSE)
# mapping$samples <- row.names(mapping)
# 
# data <- read.table(args[1], sep='\t', header=T, row.names=1,check.names=F)
# 
# #取供体的top30
# data1 <- data
# Donor <- colnames(data1)[grepl('Donor', colnames(data1))]
# # data1_t <- data1[data1$V6 != "unclassified", ]
# Donor_df <- data.frame('value'=data1[, Donor], 'kegg'=row.names(data))
# Donor_df <- Donor_df[Donor_df$kegg != "unclassified", ]
# top_s <- Donor_df$kegg[order(-Donor_df$value)[1:30]]
# 
# if (length(args) == 4){
#   data <- read.table(args[4], sep='\t', header=T, row.names=1,check.names=F)
# }
# 
# # write.table(data, "test.xls" ,sep='\t')
# #df <- data[top_s, mapping$samples[!grepl('Donor', mapping$samples)]]
# df <- data[top_s, mapping$samples]
# # write.table()
# # write.csv2(df, "top30.xls" ,sep='\t', col.names = NA)
# write.table(df, paste0(args[3], "/heatmap_top30.xls") ,sep='\t', col.names = NA)
# 
# split <- factor(mapping$Group, levels = unique(mapping$Group))
#绘图
fontsize1 = 12
fontsize2 = 10
# #行注释
# ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = NA, lwd=0),
#                                         labels_gp = gpar(color="black",
#                                                          fontsize=fontsize1,
#                                                          fontfamily="serif"),
#                                         labels = unique(split)),
#                        annotation_name_gp = gpar(col = "black",
#                                                  fontsize = 0,
#                                                  fontfamily='serif'),
#                        show_legend = c(FALSE),
#                        annotation_name_rot = 0,
#                        which = "col"
# )

heatmap_legend_param = list(direction = "vertical",
                            legend_height = unit(4, "cm"),
                            legend_position = "left",
                            title_gp=gpar(
                              fontsize=fontsize1, fontfamily='serif'),
                            labels_gp=gpar(
                              fontsize = fontsize2,fontfamily='serif')
)


#颜色
library(circlize)

col_fun = colorRamp2(c(quantile(as.matrix(df), 0.02), 
                       median(as.matrix(df)), 
                       quantile(as.matrix(df), 0.98)), 
                     c("#0B1279", "white", "#B83534"))

p <- Heatmap(df,  name=" ", col = col_fun,
             cluster_columns = FALSE, rect_gp = gpar(col = "grey", lwd = 2),
             row_names_gp = gpar(col = "black", alpha = 1, fontsize = fontsize1*0.75,
                                 fontfamily='serif'), 
             column_names_gp = gpar(fontsize = 12, fontfamily='serif'),
             column_names_rot = 0,
             column_names_centered = T,
             heatmap_legend_param = heatmap_legend_param,
             column_title_gp = gpar(fontsize = 12, fontfamily='serif'),
             column_title_rot= 20, 
             row_dend_width = unit(25, "mm")
             #column_gap = unit(1.5, "mm"),
             #column_split = split
             # top_annotation  = ha
)
pdf(paste0(args[3], "/heatmap_diff_top30_average.pdf"), width = 6.5, height = 7)
print(p)
dev.off()

