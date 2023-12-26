library(ggplot2)
library(ggsci)
library(ggnewscale)

test <- data.frame(y = runif(200,0,20),group = as.character(as.roman(1:4)),fac = LETTERS[1:5])

ggplot(test,aes(x = group,y = y))+
  geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,fill = fac),alpha = 0.008,show.legend = F)+
  #scale_fill_rickandmorty()+
  scale_fill_frontiers()+
  new_scale_fill()+#创建一个新的填充颜色比例尺,把后面的fill映射和前面的分离开,ggnewscale包
  geom_boxplot(aes(color = as.factor(group),fill = as.factor(group)),
               size = 0.8,#线条宽度
               width = 0.75,#箱体宽度
               show.legend = F)+
  geom_text(aes(x = group,y = -1,label = group))+#添加group文本标签
  geom_hline(yintercept = 0)+#添加水平线
  geom_segment(aes(x = 1.5,xend = 1.5,y = -Inf,yend = 0))+#添加分离group文本标签垂直辅助线
  geom_segment(aes(x = 2.5,xend = 2.5,y = -Inf,yend = 0))+
  geom_segment(aes(x = 3.5,xend = 3.5,y = -Inf,yend = 0))+
  facet_grid(. ~ fac,#横向排列根据fac变量分面
             switch = "x")+#默认情况下，标签显示在绘图的顶部和右侧。如果为“x”，则顶部标签将显示在底部。如果为“y”，右侧标签将显示在左侧。也可以设置为“both”。
  theme_classic()+
  scale_color_npg()+
  scale_fill_npg(alpha = 0.5)+
  labs(x = "Group",y = "")+
  theme(strip.placement = "inside",#分组标签放在面板内部
        panel.background = element_rect(color = "grey50",fill = NA),
        strip.text.y = element_blank(),#隐藏y轴分面标签
        panel.spacing = unit(0,"mm"),#间距
        axis.ticks.x = element_blank(),#隐藏x轴刻度
        axis.text.x = element_blank(),#隐藏x轴标签
        strip.background = element_rect(linewidth = 0.5))

#给每个分面fill独立颜色
pal <- c("#E64B35FF","#4d97cd","#00A087FF","#F39B7FFF","#e8c559",
         "#ea9c9d","#4DBBD5FF","#a3d393","#3C5488FF",
         "#ea9c9d","#4DBBD5FF","#a3d393","#3C5488FF",
         "#ea9c9d","#4DBBD5FF","#a3d393","#3C5488FF",
         "#ea9c9d","#4DBBD5FF","#a3d393","#3C5488FF",
         "#ea9c9d","#4DBBD5FF","#a3d393","#3C5488FF")
library(ggh4x)
ggplot(test,aes(x = group,y = y))+
  geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,fill = fac),alpha = 0.008,show.legend = F)+
  #scale_fill_rickandmorty()+
  scale_fill_frontiers()+
  new_scale_fill()+#创建一个新的填充颜色比例尺,把后面的fill映射和前面的分离开,ggnewscale包
  geom_boxplot(aes(color = as.factor(group),fill = as.factor(group)),
               size = 0.8,#线条宽度
               width = 0.75,#箱体宽度
               show.legend = F)+
  facet_nested(.~ fac + group,scales = "free_x",# 分面在一起，嵌套
               strip = strip_nested(background_x = elem_list_rect(fill = pal),by_layer_x = F))+
  theme_classic()+
  scale_color_npg()+
  scale_fill_npg(alpha = 0.5)+
  labs(x = "Group",y = "")+
  theme(strip.placement = "inside",#分组标签放在面板内部
        panel.background = element_rect(color = "grey50",fill = NA),
        strip.text.y = element_blank(),#隐藏y轴分面标签
        panel.spacing = unit(0,"mm"),#间距
        axis.ticks.x = element_blank(),#隐藏x轴刻度
        axis.text.x = element_blank(),#隐藏x轴标签
        strip.background = element_rect(linewidth = 0.5))

#定义分面间隔
library(ggh4x)
ggplot(test,aes(x = group,y = y))+
  geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,fill = fac),alpha = 0.008,show.legend = F)+
  #scale_fill_rickandmorty()+
  scale_fill_frontiers()+
  new_scale_fill()+#创建一个新的填充颜色比例尺,把后面的fill映射和前面的分离开,ggnewscale包
  geom_boxplot(aes(color = as.factor(group),fill = as.factor(group)),
               size = 0.8,#线条宽度
               width = 0.75,#箱体宽度
               show.legend = F)+
  facet_nested(.~ fac + group,scales = "free_x",space = "free_x",axes = 'x',# 分面在一起，嵌套
               strip = strip_nested(
                 by_layer_x = TRUE,bleed = TRUE,
                 background_x = list(element_part_rect(side = "b",#top (t), left (l), bottom (b) or right (r)
                                                       fill = NA, colour = "black", linewidth = 1),
                                     element_rect(colour = NA, fill = NA))))+
  theme_classic()+
  scale_color_npg()+
  scale_fill_npg(alpha = 0.5)+
  labs(x = "Group",y = "")+
  theme(strip.placement = "inside",#分组标签放在面板内部
        panel.background = element_rect(color = "grey50",fill = NA),
        strip.text.y = element_blank(),#隐藏y轴分面标签
        panel.spacing.x = unit(c(0,0,0,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0,0),"mm"),#手动设置间距
        axis.ticks.x = element_blank(),#隐藏x轴刻度
        axis.text.x = element_blank(),#隐藏x轴标签
        strip.background = element_rect(linewidth = 0.5))
