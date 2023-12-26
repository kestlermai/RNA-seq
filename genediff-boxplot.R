library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggprism)
library(ggpattern)
library(agricolae)
library(scales)
setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE140713")
load("exp.Rdata")
exp_t <- t(exp)
hubgene <- c("MX1","IFIT3","GBP1","SERPING1","IFIT1","IFI44","ISG15","IFI27")
hubgene <- c("ISG15","IFI27")
# 提取目的基因的基因表达矩阵
hubgene_exp <- exp_t[, hubgene]
exp_t_merge <- cbind(group=group_list,hubgene_exp)
exp_t_merge <- as.data.frame(exp_t_merge)
#剔除有问题样本
#exp_t_merge <- exp_t_merge[!rownames(exp_t_merge) %in% c("GSM4182382","GSM4182384","GSM4182394"),]
#exp_t_merge <- exp_t_merge[!rownames(exp_t_merge) %in% c("GSM827775"),]
# 将group列中的1和2替换为"control"和"HIV"
exp_t_merge$group <- ifelse(exp_t_merge$group == 1, "control", "PAH")
# 将group列转换为分类变量
exp_t_merge$group <- factor(exp_t_merge$group, levels = c("control", "PAH"))
#转换成长数据
melt_df <- melt(exp_t_merge,id.vars ="group", variable.name = "gene", value.name = "expression")
#多个基因箱线图
p2 <- ggplot(melt_df,aes(x=gene,y=expression),color = "black")+
  geom_boxplot_pattern(aes(fill=group), size=0.8, 
                       key_glyph=draw_key_rect, # 图例的形状
                       pattern='stripe', # 线条的样式
                       pattern_spacing = 0.01,# 线条之间的间距
                       pattern_density=0.01, # 线条的密度
                       pattern_angle= 45, 
                       position=position_dodge(width=1), # 对齐方式
                       outlier.shape = NA)+ # 异常值的形状（NA表示不显示）
  #scale_fill_nejm() + # 设置颜色映射
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))+
  stat_compare_means(aes(group=group),label="p.signif", method ="t.test",size=6)+#t.test; wilcox.test; anova; kruskal.test;
  scale_x_discrete(name="")+
  ggtitle("GSE33463")+
  labs(y = "Expression")+
  #scale_y_continuous(expand = c(0,0))+#从零开始
  theme_prism(border = T, base_rect_size = 1)+
  theme(axis.text=element_text(color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  guides(fill = guide_legend(override.aes = list(alpha = 1),
                             ncol = 1, byrow = TRUE),  # 设置图例分为4行，并显示所有图例项
         keywidth = unit(4, units = "mm"))+ # 设置图例项宽度
  theme(axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(face = "bold"),
        legend.position = "right", legend.direction = "vertical")+#horizontal,vertical
  theme(legend.position=c(0.1, 0.9))
p2
#提取nejm颜色
pal_nejm(palette = c("default"), alpha = 1)(8)
show_col(pal_nejm(palette = c("default"), alpha = 1)(8))
#单个基因箱线图+误差线
p4 <- ggplot(melt_df,aes(x=group,y=expression),color = "black")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 1),width = 0.2,
               color = "black",size = 1) +
  geom_boxplot_pattern(aes(fill=group), size=1, 
                       key_glyph=draw_key_rect, # 图例的形状
                       pattern='stripe', # 线条的样式
                       pattern_spacing = 0.01,# 线条之间的间距
                       pattern_density=0.01, # 线条的密度
                       pattern_angle= 45, 
                       position=position_dodge(width=1), # 对齐方式
                       outlier.shape = NA)+ # 异常值的形状（NA表示不显示）
  #scale_fill_aaas() + 
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))+
  stat_compare_means(aes(group=group),label="..p.signif..", method ="t.test", size=5, fontface="bold",label.x = 1.5)+
  scale_x_discrete(name="")+
  ggtitle("GSE140713")+
  labs(y = "IFI27_Expression")+
  #scale_y_continuous(expand = c(0,0))+#从零开始
  theme_prism(border = T, base_rect_size = 1)+
  theme(axis.text=element_text(color="black"))+
  guides(fill = guide_legend(override.aes = list(alpha = 1),
                             ncol = 1, byrow = TRUE),  # 设置图例分为4行，并显示所有图例项
         keywidth = unit(4, units = "mm"))+ # 设置图例项宽度
  theme(axis.title = element_text(size = 12, face = "bold"), 
        axis.text = element_text(face = "bold"),
        legend.position = "right", legend.direction = "vertical")#horizontal,vertical
p4

library(patchwork)
p <- p3 + p2 + p1 + p4
p
ggsave("C:/Users/maihuanzhuo/Desktop/diff.pdf", p, width = 30, height = 25, dpi = 1200, units = "cm")
