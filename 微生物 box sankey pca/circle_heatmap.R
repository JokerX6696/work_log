##代码1（标准化圈图）
rm(list = ls())

args <- commandArgs(T)
##安装包

# install.packages("dendextend","circlize","openxlsx")

##加载包

library(circlize)

library(dendextend)

# library(openxlsx)

setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")
# library(GlobalOptions)
# library(colorspace)
# source("R/circos.heatmap.xyf.R")
# source("R/global.R")
# source("R/plot.R")
# source("R/utils.R")
# source("R/link.R")
# source("R/zzz.R")
args<- c("KEGG.level3.abundance.xls",
         "mapping.txt",
         "circle_heatmap", "t")##样方-环境属性矩阵



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
# results1 <- res[res$pvalue < 0.05,]

if(!dir.exists(args[3])){
  dir.create(args[3])
}

write.table(res, paste0(args[3], '/KEGG_L3_circle_heatmap_diff.xls'), 
            row.names = T, sep = '\t', quote = FALSE, col.names = NA)

results1 <- res
results1[, "sum"] <- rowSums(results1[, mapping$Sample])
results1 <- results1[order(-results1$sum),]
tmp1 <- results1[results1$significant != "not significant", ]
tmp2 <- results1[results1$significant == "not significant", ][1:(50-nrow(tmp1)),]
results1 <- rbind(tmp1, tmp2)
for(i in unique(mapping$Group)){
  results1[, i] <- rowMeans(results1[, mapping$Sample[mapping$Group==i] ])
}

results1[, "ratio"] <- results1[, unique(mapping$Group)[1]]/results1[, unique(mapping$Group)[2]]
results1[results1$pvalue<0.05 & results1$ratio <1, "color"] <- "#BC3C29"
results1[results1$pvalue<0.05 & results1$ratio >1, "color"] <- "#0072B5"
results1[is.na(results1$color), "color"] <- "black"


mat1 <- results1[, c(unique(mapping$Group))]






##读入数据
# mat1 <- read.delim(args[1], row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#添加分组信息
# group <- read.delim(args[2], sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#((C:\\Users\\xyf\\Desktop\user\工作空间\\热图\\data.xlsx）替换为你的数据的路径）

# row.names(mat1)<-mat1[,1]#修改行名

# mat1<-mat1[,-1]##删除第一列，使之变为数字矩阵，绘图的数据要求为矩阵（也就是单一类型的数据矩阵，这里全为数字）

#图片标题换行
#strng_t为字符串；nrow为转为多少行；
enter_row <- function(strng_t, nrow=2){
  tmp_l <- strsplit(strng_t, " ")[[1]]
  k = 1
  step = length(tmp_l)/nrow
  res_tmp <- c()
  while(k <= length(tmp_l)){
    if(k==1){
      res_tmp <- tmp_l[k:(k+step)]
    }else{
      res_tmp <- c(res_tmp, "\n", tmp_l[k:(k+step)])
    }
    k = k + step
  }
  res_tmp <- res_tmp[!is.na(res_tmp)]
  return(paste(res_tmp, sep = " ",collapse = " "))
}

##绘图
col_fun1 = colorRamp2(c(quantile(as.matrix(mat1), 0.01), 
                        median(as.matrix(mat1)), 
                        quantile(as.matrix(mat1), 0.99)), 
                      c("#0B1279", "white", "#B83534"))

# col_fun1 = colorRamp2(c(0, 150, 300), c("blue", "white", "red"))##设置热图颜色

column_od = hclust(dist(t(mat1)))$order #对列聚类
row_od = hclust(dist(mat1))$order #对行聚类
mat1 <- mat1[row_od, column_od]

color_text <- results1[row_od, "color"]
row_names_text <- c()
for(i in row.names(mat1)){
  if (nchar(i) >= 25){
    row_names_text <- c(row_names_text, enter_row(i, 2)) #换行
  }else{
    row_names_text <- c(row_names_text, i)
  }
}
row.names(mat1) <- row_names_text


pdf(paste0(args[3], "/circle_heatmap.pdf"), width = 7, height = 7)

# circos.clear()

circos.par(gap.after = c(10))##为添加列名留出空间

circos.heatmap(mat1, ##将列聚类后重新排序的矩阵
               col = col_fun1, ##设置颜色
               cell.border = "black",
               cell.lwd = 1.2,
               bg.border = "black",
               bg.lwd = 1.2,
               rownames.col = color_text,
               rownames.font = 1,
               rownames.cex = 0.5,
               rownames.side = "outside",##行名在圈外
               cluster = TRUE,
               # cell_width = c(1,0.5),
               # rownames.family = "serif",
               dend.track.height = 0.1,
               dend.side = "inside" ##树状图在圈内
               # dend.callback = function(dend, m, si) {
               #   # when k = 1, it renders one same color for the whole dendrogram
               #   color_branches(dend, k = 4, col = 2:5)##对树状图进行着色
               # 
               # },
               # cell_width = 0.1
)
# m <- mat1
# circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.2, 
#              panel.fun = function(x, y) {
#                sector.numeric.index = CELL_META$sector.numeric.index
#                m = mat1
#                od = CELL_META$row_order
#                nr = nrow(m)
#                if (!is.null(rownames(m))) {
#                  circos.text(CELL_META$cell_middle[od], rep(0, nr), 
#                              rownames(m)[od], 
#                              cex = 0.45, 
#                              font = 1, 
#                              family = "serif",
#                              col = col_fun1[od], 
#                              facing = "clockwise", niceFacing = TRUE, 
#                              adj = c(0, 0.5))
#                }
#              })

circos.track(track.index = 2, ##将列名添加在第二个轨道（就是热图所在的环形轨道）
             
             panel.fun = function(x, y) {
               
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 
                 cn = colnames(mat1)##取得列名
                 
                 n = length(cn)
                 
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), ##x轴坐标
                             
                             1:n - convert_y(0.5, "mm"), ##y轴坐标
                             
                             labels = cn, ##输入要展示的列名
                             col = "black",
                             cex = 0.5, ##列名的大小
                             # crt = 90,
                             # srt = 90,
                             adj = c(1, 0.8),
                             family = "serif",
                             facing = "clockwise"
                 )
                 
               }
               
             }, bg.border = NA)

circos.clear()

lgd = Legend(title = "Abundance", 
             legend_height = unit(0.7, "inches"),
             title_position ="topcenter",
             labels_gp = gpar(fontsize = 10, fontfamily = "serif"),
             title_gp = gpar(fontsize = 12, fontfamily = "serif"),
             # legend_gp = gpar(fontsize = 20,fontfamily="serif"),
             col_fun = col_fun1)
grid.draw(lgd)

dev.off()
