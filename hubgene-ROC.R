library(dplyr)
library(pROC)
library(ggplot2)
library(ggprism)
library(ggsci)
library(verification)
library(RColorBrewer)
setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE703")
load("exp.Rdata")
exp_t <- t(exp)
exp_t_merge <- cbind(group=group_list,exp_t)
exp_t_merge <- as.data.frame(exp_t_merge)
#剔除有问题样本
#exp_t_merge <- exp_t_merge[!rownames(exp_t_merge) %in% c("GSM827775"),]
# 将group列中的1和2替换为"control"和"HIV"
exp_t_merge$group <- ifelse(exp_t_merge$group == 1, "control", "PAH")
# 将group列转换为分类变量
exp_t_merge$group <- factor(exp_t_merge$group, levels = c("control", "PAH"))
# 提取hub基因
hubgene <- c("MX1","IFIT3","GBP1","SERPING1","IFIT1","IFI44","ISG15","IFI27")
hubgene <- c("ISG15","IFI27")
#hubgene <- c("MX1","SERPING1","IFIT1","ISG15","IFI27")
#hubgene <- "IFI27"

# 创建一个空数据框，用于存储结果
result <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(result) <- c("gene", "AUC", "P.value","95%CI")
# 循环计算每列对应的AUC值和P值
for(gene in hubgene) {
  roc <- roc(group ~ exp_t_merge[[gene]], exp_t_merge,levels = c("control", "PAH"))
  auc <- roc$auc
  ci <- ci(roc, level = 0.95)  # add this line to compute confidence interval
  row <- data.frame(gene = gene, AUC = auc, Lower_CI = ci[1], Upper_CI = ci[2])  # add Lower_CI and Upper_CI columns
  if(length(levels(exp_t_merge$group)) == 2){
    exp_t_merge$result <- ifelse(exp_t_merge$group == "control", 0, 1)
    pval <- roc.area(exp_t_merge$result, roc$predictor)$p.value
    row$`P.value` <- pval
  }
  result <- rbind(result, row)
}
write.csv(result, "ROC-result.csv", row.names = FALSE)  

# 循环计算每个基因的ROC曲线和95%CI
roc_list <- list()
auc_list <- list()
CI_list <- list()
for (gene in hubgene) {
  roc <- roc(group ~ exp_t_merge[[gene]], exp_t_merge, levels = c("control", "PAH"))
  auc <- roc$auc
  ci <- ci(roc, level = 0.95) # 计算95%CI
  exp_t_merge$result <- ifelse(exp_t_merge$group == "control", 0, 1)
  pval <- roc.area(exp_t_merge$result, roc$predictor)$p.value
  roc_list[[gene]] <- roc
  auc_list[[gene]] <- auc
  CI_list[[gene]] <- ci
}

labels <- paste0(names(roc_list), ", AUC = ", format(round(sapply(auc_list, function(x) round(x, 3)), 3), nsmall = 3),
                 ", 95%CI (", format(round(sapply(CI_list, function(x) round(x[1], 3)), 3), nsmall = 3), "-",
                 format(round(sapply(CI_list, function(x) round(x[2], 3)), 3), nsmall = 3), ")")
head(labels)

#绘图
p1 <- ggroc(roc_list, legacy.axes = TRUE,linewidth = 1)+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), #添加对角线
               color="darkgrey",linewidth = 1,linetype="dashed")+
  ggtitle("GSE703") + xlab("1 - Specificity") + ylab("Sensitivity") +
  #theme_bw()+
  theme_prism(border = T, base_rect_size = 1) +
  #theme_minimal()+
  scale_color_brewer(labels = labels, type = "qual", palette = "Set1")+#RColorBrewer::display.brewer.all() 
  theme(legend.position = c(0.70,0.25),legend.key.width = unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol = 1, keyheight = unit(0.8, "lines"), 
                              label.position = "right",label.hjust = 0))+
  theme(plot.title = element_text(face = "bold",hjust = 0.5), 
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(face = "bold", size = 10), legend.title = element_blank(),
        legend.text = element_text(face = "bold",size = 10))
p1
library(patchwork)
p <- p4 + p2 + p3 + p1
p
ggsave("C:/Users/maihuanzhuo/Desktop/roc.pdf", p, width = 30, height = 25, dpi = 1200, units = "cm")
