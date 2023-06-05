#!/usr/bin/env Rscript

rm(list = ls())
# library(docopt)
# setwd("D:/02.项目/二代宏基因组/DOE20225714_邓哲_售后_20230523")


suppressPackageStartupMessages(library("optparse"))

# "Usage: alpha_estimator_boxplot.r  -i <file> --output <dir> -g <file> [--grouporder <string>]
# Options:
#    -i , --input <file>            alpha_estimator_summary.xls
#    -g , --groupfile <file>        group info
#    -o , --output <dir>            output dir
#    --grouporder <file>            give the order of group, \'group1,group2,group3\' [default: NA]" -> doc

# opts                     <- docopt(doc, version = 'alpha 多样性指数 箱体图\n')
# input                   <- opts$input
# groupfile                <- opts$groupfile
# grouporder               <- opts$grouporder
# output                   <- opts$output

option_list = list(
make_option(c("-i", "--input"), type = "character", default = NULL,
help = "input file name", metavar = "character"),
make_option(c("-m", "--mapping"), type = "character", default = NULL,
help = "mapping file name", metavar = "character"));
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


# input = "/home/lihongjie/pipeline/16S/2021_12_6/alpha_analysis/alpha_index_min/alpha_estimator_summary.xls"
# groupfile = "/home/lihongjie/pipeline/16S/2021_12_6/alpha_analysis/alpha_index_min/sample.group.xls"
# output = "/home/lihongjie/pipeline/16S/2021_12_6/alpha_analysis/alpha_index_min/"

input = opt$input
mapping = opt$mapping
# input = "alpha_estimator_summary.xls"
# mapping = "mapping.txt"
output = "./"
grouporder = "NA"

library("ggplot2")
library("dplyr")
# library("ggpubr")
library("ggsignif")

data = read.table(input, sep="\t", header=T, check.name = FALSE, quote = "")
group = read.table(mapping, sep="\t", header=T, check.name = FALSE, quote = "", comment.char ="")
# group = group[, c(1, 4)]
colnames(group) = c("Samples","Group")

if (grouporder != 'NA'){
	custom_level <- unlist(strsplit(grouporder,split=","))
}else{
	custom_level <- unique(group$Group)
}

out = merge(x=data, y=group, by = "Samples", all = TRUE)
alpha_index = colnames(out)

group_num = length(unique(as.character((group$Group))))
if (group_num < 3){
	width = 1.4 + 1.0*group_num
}else{
	width = 1.4 + 0.6*group_num
}
height = 3.4 + 0.2*group_num

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
taxon = alpha_index

for (i in 2:(length(alpha_index)-1)){
	plotdata = data.frame(Group=out$Group, abundance=out[,i])
	plotdata$Group<-factor(plotdata$Group,levels=custom_level)
	tongji <- plotdata %>% group_by(Group) %>% summarise(value = sum(abundance))
	if (length(tongji$value[tongji$value %in% 0]) > 0) next # 一整组数值为零，跳过

	res <- pairwise.t.test(plotdata$abundance, plotdata$Group, p.adjust.method='none')
	# pValue <- c(res$p.value)[!is.na(c(res$p.value))]
	# my_annotations <- as.character(symnum(pValue, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***",'**','*','')))
	# my_annotations1 <- my_annotations[my_annotations != ""]

	pValue <- c(res$p.value)[!is.na(c(res$p.value))]
	my_annotations <- as.character(symnum(pValue, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***",'**','*','')))
	pValue_bak <- pValue[my_annotations != ""]
	my_annotations1 <- my_annotations[my_annotations != ""]
	
	
	wtq <- levels(plotdata$Group)
	lis <- combn(wtq,2)
	my_comparisons <- tapply(lis, rep(1:ncol(lis), each=nrow(lis)), function(i) i)
	my_comparisons1 <- my_comparisons[my_annotations != ""]

	my_levels <- custom_level


	p <- ggplot(data=plotdata, aes(x=Group,y=abundance,fill=Group))+
	  stat_boxplot(geom = "errorbar",width=0.2, alpha = 1, linetype = 1, lwd = 1) + 
	  geom_boxplot(
	    colour='black',
	    alpha = 1,
	    linetype = 1,
	    size = 0.8,
	    width = 0.7,
	    lwd = 1,
	    notch = T,
	    notchwidth = 0.7,
	    outlier.alpha = 0)+
	  # geom_jitter(
	  #   aes(shape=Group,colour=Group),
	  #   alpha = 1, show.legend = FALSE,
	  #   size = 1.5)+
	  geom_jitter(aes(color=Group), shape = 21, size=3, alpha = 0.5, 
	              colour='black',
	              position = position_jitterdodge(jitter.width = 0.85,
	                                              dodge.width = 0.75), #width =0.25,
	              show.legend = F)+ ##添加散点
	  scale_fill_manual(limits=my_levels, 
	                     values=mycol2[1:length(my_levels)])+
	  scale_color_manual(limits=my_levels, 
	                    values=mycol2[1:length(my_levels)])+
	  theme_classic(base_line_size = 1)+
	  coord_cartesian(clip = 'off')+
	  theme(
	    text=element_text(family="serif"),  #设置字体
	    # panel.background = element_rect(fill=NULL),
	    panel.border = element_rect(size=1, colour='black', fill=NA),
	    plot.title = element_text(size = 12,
	                              colour = "black",
	                              face = "bold",
	                              hjust = 0.5),
	    axis.title.y = element_text(size = 13, 
	                                color = "black",
	                                face = "bold", 
	                                vjust = 1.9, 
	                                hjust = 0.5, 
	                                angle = 90),
	    legend.title = element_text(color="black", # 修改图例的标题
	                                size=13, 
	                                face="bold"),
	    legend.text = element_text(color="black", # 设置图例标签文字
	                               size = 11),
	    axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
	                               color = "black", # 颜色
	                               vjust = 1, # 位置
	                               hjust = 0.5, 
	                               angle = 0), #角度
	    axis.text.y = element_text(size = 11,  
	                               color = "black",
	                               vjust = 0.5, 
	                               hjust = 0.5, 
	                               angle = 0),
	    axis.ticks.x= element_line(size=1),
	    axis.ticks.y= element_line(size=1), 
	    axis.line.x = element_line(size = 1), 
	    axis.line.y = element_line(size = 1))
	
	#添加P值
	if (length(my_annotations1) > 0){
	  p = p + geom_signif(annotations = '',
	                      comparisons = my_comparisons1,
	                      step_increase = 0.1,
	                      vjust=.2, colour='gray20', 
	                      tip_length=0.015)
	  
	  library(grid)
	  max_v = max(plotdata$abundance)
	  min_v = min(plotdata$abundance)
	  scale_v = (max_v-min_v)/100
	  
	  for (k in 1:length(my_comparisons1)){
	    x_t = which(my_levels %in% my_comparisons1[[k]])
	    if (taxon[i] == 'goods_coverage'){
	      y_p = max_v+k*scale_v*7+scale_v*1
	    }else{
	      y_p = max_v*0.9759999+k*max_v*0.09400003-scale_v*1
	    }
	    
	    p <- p + annotation_custom(
	      grob = textGrob(label = substitute(italic(p)~'='~x, 
	                                         list(x=round(pValue_bak[k], 4))), 
	                      just='center',
	                      hjust='center',
	                      vjust='bottom',
	                      gp = gpar(
	                        cex = 0.7, 
	                        fontfamily="serif")),
	      ymin = y_p,      # a=0.0662094,max_v=0.7043551,ymin=0.7536599
	      ymax = y_p,
	      xmin = x_t[1],         # Note: The grobs are positioned outside the plot area
	      xmax = x_t[2])
	  }
	  p <- p + scale_y_continuous(limits = c(min_v-scale_v*1,y_p+scale_v*3))
	}
	
	#图片宽度
	modify_name <- function(x, cex = 1.043) {
	  # par()
	  x <- as.character(x)
	  w <- strwidth(x, family = "serif", units = "inches", cex = cex)
	  return(w)
	}
	pdf(NULL)
	lwidth <- max(vapply(unique(plotdata$Group), modify_name, 0, USE.NAMES = F))
	dev.off()
	group_num = length(levels(plotdata$Group))
	if (group_num < 3){
	  width = 2.7 + 1.0*group_num+lwidth
	}else{
	  width = 1.18 + 0.6*group_num+lwidth
	}
	#图片高度
	pdf(NULL)
	lheight <- max(vapply(unique(plotdata$Group), modify_name, 0, USE.NAMES = F))
	dev.off()
	height = 3.9 + lheight*sqrt(2)/2
	#图片标题
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
	p <- p+labs(x="",y=substitute(x, list(x=taxon[i])))
	# #添加标题，过长换行
	# pdf(NULL)
	# taxon_name <- modify_name(taxon[i], cex = 1.043)
	# dev.off()
	# if (taxon_name > width*0.85){
	#   p <- p+labs(title=enter_row(taxon[i], ceiling(taxon_name/(width*0.85))),x="",y="Abundance")
	# }else{
	#   # library(grDevices)
	#   p <- p+labs(x="",y=substitute(paste(italic(x), " (RA)"), list(x=taxon[i])))
	# }
	# 
	
  #if (length(my_annotations1) > 0) p = p + geom_signif(annotations = my_annotations1, comparisons = my_comparisons1, step_increase = 0.1, vjust=.2, colour='gray20', tip_length=0.015)
	pdf(file=paste0(output, "/alpha_estimator_boxplot", "_", alpha_index[i], ".pdf"), width=width, height=height)
	print(p)
	dev.off()
}
# dev.off()
