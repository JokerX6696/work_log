rm(list=ls())
setwd('D:\\desk\\DZOE2023041442-王石-数据高级分析')
library("VennDiagram")
library(grid)
library(futile.logger)
file <- c('F2-1_0.8_stat.xls','F2-2_0.8_stat.xls')
df1 <- read.table(file[1],sep = '\t',header = F,row.names = 1)
df2 <- read.table(file[2],sep = '\t',header = F,row.names = 1)
#df3 <- read.table(file[3],sep = '\t',header = F,row.names = 1)
df = df1 + df2# +df3
names(df) <- 'reads_num'
N12 <- df$reads_num[1] + df$reads_num[7]
AREA1 <- df$reads_num[3] + N12
AREA2 <- df$reads_num[4] + N12 

out_file <- 'F2_venn_group.pdf'
pdf(file = out_file,width = 6,height = 6)

venn.plot <- draw.pairwise.venn(
  area1 = AREA1,
  area2 = AREA2,
  scaled = FALSE,
  cross.area = N12,
  category = c("Cyprinus carpio", "Megalobrama amblycephala"),
  fill = c("#00468BFF", "#42B540FF"),
  cat.pos = c(180,0),
  lty = "blank",
  cex = 1,
  cat.cex = 1.5,
  cat.col = c("#00468BFF", "#42B540FF"),
  margin = 0.15,
  cex.label = 0.5
  #label.col = "black"
);
#grid.draw(venn.plot);#画图展示
dev.off()
png(file = sub("pdf","png",out_file),width = 600,height = 600)
#par(mar=c(3,3,2,3),oma=c(1,1,1,1))
venn.plot <- draw.pairwise.venn(
  area1 = AREA1,
  area2 = AREA2,
  scaled = FALSE,
  cross.area = N12,
  category = c("Cyprinus carpio", "Megalobrama amblycephala"),
  fill = c("#00468BFF", "#42B540FF"),
  cat.pos = c(180,0),
  lty = "blank",
  cex = 1,
  cat.cex = 1.5,
  cat.col = c("#00468BFF", "#42B540FF"),
  margin = 0.15,
  cex.label = 0.5
  #label.col = "black"
);
#grid.draw(venn.plot);#画图展示
dev.off()
