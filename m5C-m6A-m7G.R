library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggprism)
library(ggpattern)
library(reshape2)
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE140713")
load("exp.Rdata")

###构建m6A基因集
m6A_genes <- c("METTL3","METTL14","METTL16","RBM15","RBM15B","VIRMA","ZC3H13","FTO","ALKBH3","ALKBH5",
               "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1","IGF2BP2","IGF2BP3","HNRNPA2B1","HNRNPC")

#从表达矩阵中筛选出m6A相关基因的表达量
m6A_exp <- exp[rownames(exp) %in% m6A_genes, ]
m6A_exp <- as.data.frame(m6A_exp)

#合并分组信息
m6A_exp <- t(m6A_exp)
m6A_exp_merge <- cbind(group = group_list, m6A_exp)
m6A_exp_merge <- as.data.frame(m6A_exp_merge)
# 将group列中的1和2替换为"control"和"HIV"
m6A_exp_merge$group <- ifelse(m6A_exp_merge$group == 1, "control", "PAH")
# 将group列转换为分类变量
m6A_exp_merge$group <- factor(m6A_exp_merge$group, levels = c("control", "PAH"))
#将矩阵数据转换成长数据进行可视化
melt_df <- melt(m6A_exp_merge, id.vars ="group", variable.name = "gene", value.name = "expression")

###构建m7G基因集
m7G_genes <- c("AGO2","CYFIP1","DCP2","DCPS","EIF3D","EIF4A1","EIF4E","EIF4E2","EIF4E3","EIF4G3",
               "GEMIN5","IFIT5","LARP1","LSM1","METTL1","NCBP1","NCBP2","NCBP2L","NCBP3","NSUN2","NUDT10",
               "NUDT11","NUDT16","NUDT3","NUDT4","SNUPN","WDR4")
#从表达矩阵中筛选出m6A相关基因的表达量
m7G_exp <- exp[rownames(exp) %in% m7G_genes, ]
m7G_exp <- as.data.frame(m7G_exp)
#合并分组信息
m7G_exp <- t(m7G_exp)
m7G_exp_merge <- cbind(group = group_list, m7G_exp)
m7G_exp_merge <- as.data.frame(m7G_exp_merge)
#剔除有问题样本
#m7G_exp_merge <- m7G_exp_merge[!rownames(m7G_exp_merge) %in% c("GSM827775"),]
# 将group列中的1和2替换为"control"和"HIV"
m7G_exp_merge$group <- ifelse(m7G_exp_merge$group == 1, "control", "HIV")
# 将group列转换为分类变量
m7G_exp_merge$group <- factor(m7G_exp_merge$group, levels = c("control", "HIV"))
#将矩阵数据转换成长数据进行可视化
melt_df <- melt(m7G_exp_merge, id.vars ="group", variable.name = "gene", value.name = "expression")

###构建m5C基因集
m5C_genes <- c("NOP2","NSUN2","NSUN3","NSUN4","NSUN5","NSUN6","NSUN7","DNMT1","DNMT2","DNMT3A","DNMT3B",
               "TRDMT1","TET2","ALYREF")
#从表达矩阵中筛选出m6A相关基因的表达量
m5C_exp <- exp[rownames(exp) %in% m5C_genes, ]
m5C_exp <- as.data.frame(m5C_exp)
#合并分组信息
m5C_exp <- t(m5C_exp)
m5C_exp_merge <- cbind(group = group_list, m5C_exp)
m5C_exp_merge <- as.data.frame(m5C_exp_merge)
#剔除有问题样本
#m5C_exp_merge <- m5C_exp_merge[!rownames(m5C_exp_merge) %in% c("GSM827775"),]
# 将group列中的1和2替换为"control"和"HIV"
m5C_exp_merge$group <- ifelse(m5C_exp_merge$group == 1, "control", "PAH")
# 将group列转换为分类变量
m5C_exp_merge$group <- factor(m5C_exp_merge$group, levels = c("control", "PAH"))
#将矩阵数据转换成长数据进行可视化
melt_df <- melt(m5C_exp_merge, id.vars ="group", variable.name = "gene", value.name = "expression")

p2 <- ggplot(melt_df,aes(x=gene,y=expression),color = "black")+
  geom_boxplot_pattern(aes(fill=group), size=0.8, 
                       key_glyph=draw_key_rect, # 图例的形状
                       pattern='stripe', # 线条的样式
                       pattern_spacing = 0.01,# 线条之间的间距
                       pattern_density=0.01, # 线条的密度
                       pattern_angle= 45, 
                       position=position_dodge(width=1), # 对齐方式
                       outlier.shape = NA)+ # 异常值的形状（NA表示不显示）
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))+ # 设置颜色映射
  stat_compare_means(aes(group=group),label="p.signif", method ="t.test",size=5)+#t.test; wilcox.test; anova; kruskal.test;
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
        legend.position = "right", legend.direction = "vertical")#horizontal,vertical
p2

library(patchwork)
p <- p1/p2
p
ggsave("C:/Users/11146/Desktop/m7G.pdf", p, width = 30, height = 30, dpi = 1200, units = "cm")

#与IFI27相关性分析
#添加IFI27加入表达矩阵中
IFI27_exp <- exp[rownames(exp) %in% "IFI27", ]
IFI27_exp <- as.data.frame(IFI27_exp)
colnames(IFI27_exp) <- "IFI27"
IFI27_exp <- t(IFI27_exp)
m5C_exp <- t(m5C_exp)
m5C_exp_IFI27 <- rbind(m5C_exp, IFI27_exp)

y <- as.numeric(m5C_exp_IFI27["IFI27", ])
#创建一个空的数据框用于存储结果
result_df <- data.frame(Gene = character(), Correlation = numeric(), pValue = numeric(), stringsAsFactors = F)
#循环基因集
for (i in m5C_genes) {##rownames(m5C_exp_IFI27)
  x <- as.numeric(m5C_exp_IFI27[i, ])
  df <- data.frame(x, y)
  cor.test_result <- cor.test(x, y, method = "pearson")#pearson，spearman
  cor_value <- cor.test_result$estimate
  p_value <- cor.test_result$p.value
  if (cor_value > 0.3 && p_value < 0.05) {
    result_df <- rbind(result_df, data.frame(Gene = i, Correlation = cor_value, pValue = p_value, stringsAsFactors = F))
  }
}
print(result_df$Gene)


library(ggpubr)
#install.packages("ggExtra")
library(ggExtra)
###m7G共享有DCPS和NCBP1
x <- as.numeric(m5C_exp_IFI27["LARP1",])
y <- as.numeric(m5C_exp_IFI27["IFI27",])
df <- as.data.frame(cbind(x, y))
cor.test <- cor.test(x, y, method = "pearson")#pearson，spearman
cor <- cor.test$estimate
pValue <- cor.test$p.value
#绘图
p1 <- ggplot(df, aes(x, y)) + 
  labs(x = "LARP1",y = "IFI27")+
  geom_point()+ 
  geom_smooth(method = "lm", formula = y ~ x) + 
  theme_bw()+
  theme(axis.title = element_text(size = 12))+
  stat_cor(method = 'pearson', aes(x = x, y = y))
p2 <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
ggsave("C:/Users/11146/Desktop/hiv_LARP1.pdf", p2, width = 15, height = 12, dpi = 1200, units = "cm")
