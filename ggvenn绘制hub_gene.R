library(ggvenn)
library(tidyverse)
library(ggtext)
library(magrittr)
library(ggpubr)
library(cowplot)

setwd("C:/Users/maihuanzhuo/Desktop")

#导入数据
A <- read_tsv("A_diff.txt",col_names = F)#数据格式：一串gene name，没有列名
B <- read_tsv("B_diff.txt",col_names = F)
C <- read_tsv("WGCNA.txt",col_names = F)

#绘制venn图
p1 <- list(A = A$X1,B = B$X1,C = C$X1) %>% 
  ggvenn(digits = 1,#保留小数点后位数
         stroke_color = "white",#圆边颜色
         fill_color = c("#1E90FF", "#FF8C00","#4DAF4A"),
         set_name_color = c("#1E90FF","#FF8C00","#4DAF4A"),
         show_percentage = T,#展示百分比
         show_elements = F,#显示集合内元素
         label_sep = ',')#元素分隔符
p1

xx <- list(A = A$X1,B = B$X1,C = C$X1)
# 使用Reduce函数来计算列表中所有向量的交集
Intersect <- function(x) {
  Reduce(intersect, x)
}

# 使用Reduce函数来计算列表中所有向量的并集
Union <- function(x) {
  Reduce(union, x)
}

# 计算两个列表的差集
diff <- function(x, y) {
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

#绘制hub-gene
p2 <- Intersect(xx) %>% as.data.frame() %>% set_colnames("hub-gene") %>% 
  ggtexttable(rows = NULL, theme = ttheme("lBlueWhite"))
p2

#拼接
p1 %>% ggdraw() + draw_plot(p2,scale = 0.008,x = 0.62,y = 0.27,width = 0.5,height = 0.1)
ggsave(filename = "hub-gene-venn.pdf", width = 20, height = 14, units = "cm", dpi = 600)
