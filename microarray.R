library(GEOquery)
library(geneExpressionFromGEO)
library(limma)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(cowplot)

setwd("C:/Users/maihuanzhuo/Desktop/RNA-seq/microarray")

# GSE110224
## 下载数据
eSetGSE110224 <- getGEO("GSE110224", destdir = '.', getGPL = F)
# (1)提取表达矩阵exp
exp_1 <- exprs(eSetGSE110224[[1]])
# (2)提取临床信息
pd_1 <- pData(eSetGSE110224[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd_1),colnames(exp_1));p
if(!p) exp1 = exp1[,match(rownames(pd_1),colnames(exp_1))]
# (4)提取芯片平台编号
gpl <- eSetGSE110224[[1]]@annotation

## 创建分组信息
group_list_1 <- factor(ifelse(str_detect(pd_1$title, "normal"), "normal", "tumor"),levels = c("normal", "tumor"))

# 转换探针，提取gene symbol
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL570')# 没找到
# 下载GEO上的探针文件
GPL570 <- getGEO('GPL570', destdir = ".")
GPL570_anno <- Table(GPL570)
colnames(GPL570_anno)
# 提取对应gene symbol
anno <- GPL570_anno[, c("ID", "Gene Symbol")]
colnames(anno)[colnames(anno) == "Gene Symbol"] <- "symbol"
# 提取列中的第一个基因名称
anno$symbol <- sapply(strsplit(anno$symbol, "///"), "[", 1)
### 匹配gene id
exprSet_1 <- exp_1 %>% 
  data.frame() %>% 
  mutate(ID = rownames(exp_1)) %>% 
  inner_join(anno, by = "ID") %>%
  select(-ID) %>% #去掉多余信息
  select(symbol, everything()) %>% #重新排列，
  mutate(rowMean = rowMeans(. [grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  filter(symbol != "---") %>%
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol, .keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")

# 绘制质量箱式图
group_list_char_1 <- as.character(group_list_1)# 转换为字符向量
color_map_1 <- c(tumor = "#EF7C2C", normal = "#2F56A6")# 定义颜色映射
group_colors_1 <- color_map_1[group_list_char_1]# 生成颜色向量
# 
pdf("boxplot-GSE110224.pdf", width = 8, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE110224 (CRC)")
boxplot(exprSet_1, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 2,
        col = group_colors,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
dev.off()
## 很明显不齐，得标化一下
exprSet_1_norm <- normalizeBetweenArrays(exprSet_1)
# 
pdf("boxplot-GSE110224.pdf", width = 8, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste("GSE110224")
boxplot(exprSet_1_norm, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 1.5,
        col = group_colors,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("tumor","normal"), fill = group_colors, border = "black", bg = "white", cex = 0.8)
dev.off()

## GSE81558
## 下载数据
eSetGSE81558 <- getGEO("GSE81558", destdir = '.', getGPL = F)
# (1)提取表达矩阵exp
exp_2 <- exprs(eSetGSE81558[[1]])
# (2)提取临床信息
pd_2 <- pData(eSetGSE81558[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd_2),colnames(exp_2));p
if(!p) exp_2 = exp_2[,match(rownames(pd_2),colnames(exp_2))]
# (4)提取芯片平台编号
gpl <- eSetGSE81558[[1]]@annotation
##
pd_2_clean <- pd_2 %>% 
  filter(`disease state:ch1` == "normal"| `disease state:ch1` == "Primary Colorectal Tumor")
## 创建分组信息
group_list_2 <- factor(ifelse(str_detect(pd_2_clean$`disease state:ch1`, "normal"), "normal","tumor"), 
                       levels = c("normal", "tumor"))
# 转换探针，提取gene symbol
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL15207')# 没找到
# 下载GEO上的探针文件
GPL15207 <- getGEO('GPL15207', destdir = ".")
GPL15207_anno <- Table(GPL15207)
colnames(GPL15207_anno)
# 提取对应gene symbol
anno <- GPL15207_anno[, c("ID", "Gene Symbol")]
colnames(anno)[colnames(anno) == "Gene Symbol"] <- "symbol"
# 提取列中的第一个基因名称
anno$symbol <- sapply(strsplit(anno$symbol, "///"), "[", 1)
### 匹配gene id
exprSet_2 <- exp_2 %>% 
  data.frame() %>% 
  mutate(ID = rownames(exp_2)) %>% 
  inner_join(anno, by = "ID") %>%
  select(-ID) %>% #去掉多余信息
  select(symbol, everything()) %>% #重新排列，
  mutate(rowMean = rowMeans(. [grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  filter(symbol != "---") %>%
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol, .keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
## 挑选非转移样本
exprSet_2 <- exprSet_2[, rownames(pd_2_clean)]
exprSet_2_norm <- normalizeBetweenArrays(exprSet_2)
# 绘制质量箱式图
group_list_char_2 <- as.character(group_list_2)# 转换为字符向量
color_map_2 <- c(tumor = "#EF7C2C", normal = "#2F56A6")# 定义颜色映射
group_colors_2 <- color_map_2[group_list_char_2]# 生成颜色向量
# 
pdf("boxplot-GSE81558.pdf", width = 8, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste("GSE81558")
boxplot(exprSet_2_norm, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 1.5,
        col = group_colors,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("tumor","normal"), fill = group_colors, border = "black", bg = "white", cex = 0.8)
dev.off()


## GSE32323
## 下载数据
eSetGSE32323 <- getGEO("GSE32323", destdir = '.', getGPL = F)
# (1)提取表达矩阵exp
exp_3 <- exprs(eSetGSE32323[[1]])
# (2)提取临床信息
pd_3 <- pData(eSetGSE32323[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd_3), colnames(exp_3));p
if(!p) exp_3 = exp_3[,match(rownames(pd_3), colnames(exp_3))]
# (4)提取芯片平台编号
gpl <- eSetGSE32323[[1]]@annotation
##
pd_3_clean <- pd_3 %>% 
  filter(source_name_ch1 == "surgically resected material")
## 创建分组信息
group_list_3 <- factor(ifelse(str_detect(pd_3_clean$title, "normal"), "normal","tumor"), 
                       levels = c("normal", "tumor"))
# 转换探针，提取gene symbol
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL570')# 没找到
# 下载GEO上的探针文件
GPL570 <- getGEO('GPL570', destdir = ".")
GPL570_anno <- Table(GPL570)
colnames(GPL570_anno)
# 提取对应gene symbol
anno <- GPL570_anno[, c("ID", "Gene Symbol")]
colnames(anno)[colnames(anno) == "Gene Symbol"] <- "symbol"
# 提取列中的第一个基因名称
anno$symbol <- sapply(strsplit(anno$symbol, "///"), "[", 1)
### 匹配gene id
exprSet_3 <- exp_3 %>% 
  data.frame() %>% 
  mutate(ID = rownames(exp_3)) %>% 
  inner_join(anno, by = "ID") %>%
  select(-ID) %>% #去掉多余信息
  select(symbol, everything()) %>% #重新排列，
  mutate(rowMean = rowMeans(. [grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  filter(symbol != "---") %>%
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol, .keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
## 挑选非转移样本
exprSet_3 <- exprSet_3[, rownames(pd_3_clean)]
exprSet_3_norm <- normalizeBetweenArrays(exprSet_3)
# 绘制质量箱式图
group_list_char_3 <- as.character(group_list_3)# 转换为字符向量
color_map_3 <- c(tumor = "#EF7C2C", normal = "#2F56A6")# 定义颜色映射
group_colors_3 <- color_map_3[group_list_char_3]# 生成颜色向量
# 
pdf("boxplot-GSE32323.pdf", width = 8, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste("GSE32323")
boxplot(exprSet_3_norm, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 1.5,
        col = group_colors,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("tumor","normal"), fill = group_colors, border = "black", bg = "white", cex = 0.8)
dev.off()

##### 合并矩阵，microarray只有平台相同以及探针注释文件相同才能合并，所以这里考虑合并1和3
# 不然，则只能通过取交集进行RRA分析
exp_merge <- cbind(exprSet_1_norm, exprSet_3_norm)
dim(exp_merge)# 23348    68
group_list_merge <- c(group_list_1, group_list_3)
table(group_list_merge)
# normal  tumor 
# 33     35 
table(group_list_1) # 34
table(group_list_3) # 34
batchType <- c(rep("GSE110224", 34), rep("GSE32323", 34))
color_map_4 <- c("GSE110224" = "#974F9E", "GSE32323" = "#50AE4A")# 定义颜色映射
group_colors_4 <- color_map_4[batchType]# 生成颜色向量
## 
pdf("boxplot-merge_before.pdf", width = 10, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste("Original")
boxplot(exp_merge, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 1.5,
        col = group_colors_4,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("GSE110224","GSE32323"), # 有毒为啥legend匹配不上颜色，后面改fill就可以了神奇
       fill = c("#974F9E","#50AE4A"), border = "black", bg = "white", cex = 0.8)
dev.off()

### PCA
batchType <- factor(batchType, levels = c("GSE110224", "GSE32323"))
library(FactoMineR)
library(factoextra)
dat.pca <- PCA(t(exp_merge), graph = F)
pca_plot_before <- fviz_pca_ind(dat.pca,
                                geom.ind = "point", # show points only (nbut not "text")
                                pointshape = 21,
                                col.ind = batchType, 
                                fill.ind = batchType,
                                palette = c("#974F9E","#50AE4A"),
                                addEllipses = T, # Concentration ellipses
                                legend.title = "Group",
                                title = "Before Correction")+
  # theme_bw()+
  ggprism::theme_prism(border = T) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14),
        legend.background = ggfun::element_roundrect(color = "black", linetype = 1, size = 0.8),
        plot.title = element_text(hjust = 0.50, size = 15),
        # panel.border = element_rect(color = 'black', linewidth = 1.0),
        axis.ticks = element_line(color = 'black', linewidth = 0.8))
pca_plot_before

### 矫正可以用combat也可用limma的removeBatchEffect()
exp_merge_after <- removeBatchEffect(exp_merge, batch = batchType, group = group_list_merge)
# combat去批次
exp_merge_after <- sva::ComBat(dat = as.matrix(exp_merge), batch = batchType, mod = model.matrix( ~ group_list_merge))
# 感觉combat去批次更好
dat.pca <- PCA(t(exp_merge_after), graph = F)
pca_plot_after <- fviz_pca_ind(dat.pca,
                               geom.ind = "point", # show points only (nbut not "text")
                               pointshape = 21,
                               col.ind = batchType, 
                               fill.ind = batchType,
                               palette = c("#974F9E","#50AE4A"),
                               addEllipses = T, # Concentration ellipses
                               legend.title = "Group",
                               title = "After Correction")+
  # theme_bw()+
  ggprism::theme_prism(border = T) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14),
        legend.background = ggfun::element_roundrect(color = "black", linetype = 1, size = 0.8),
        plot.title = element_text(hjust = 0.50, size = 15),
        # panel.border = element_rect(color = 'black', linewidth = 1.0),
        axis.ticks = element_line(color = 'black', linewidth = 0.8))
pca_plot_after
p_merge <- ggpubr::ggarrange(pca_plot_before, pca_plot_after, common.legend = T, legend = "bottom")
p_merge
ggsave("./PCA_batch_merge.pdf", p_merge, width = 13, height = 6, dpi = 600)

### 再补一个箱式图
pdf("boxplot-merge_after.pdf", width = 10, height = 5)
par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
title <- paste("Batch Corrected") # original ， Batch Corrected
boxplot(exp_merge_after, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = title, cex.main = 1.5,
        col = group_colors_4,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("GSE110224","GSE32323"), # 有毒为啥legend匹配不上颜色，后面改fill就可以了神奇
       fill = c("#974F9E","#50AE4A"), border = "black", bg = "white", cex = 0.8)
dev.off()

### 合并箱式图
pdf("boxplot-merge.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
# par(mar = c(7, 4, 2, 1))#底部、左侧、顶部和右侧的边距大小
boxplot(exp_merge, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = "Original", cex.main = 1.5,
        cex.axis = 0.5, col = group_colors_4,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("GSE110224","GSE32323"), # 有毒为啥legend匹配不上颜色，后面改fill就可以了神奇
       fill = c("#974F9E","#50AE4A"), border = "black", bg = "white", cex = 0.8)
boxplot(exp_merge_after, boxwex = 0.9,#箱子宽度
        notch = T, # notch表示是否在箱子上添加缺口
        main = "Batch Corrected", cex.main = 1.5,
        cex.axis = 0.5, col = group_colors_4,
        outline = FALSE, # outline确定是否显示异常值
        las = 2) # las设置了轴标签的方向（2表示垂直）
legend("topright", legend = c("GSE110224","GSE32323"), # 有毒为啥legend匹配不上颜色，后面改fill就可以了神奇
       fill = c("#974F9E","#50AE4A"), border = "black", bg = "white", cex = 0.8)
dev.off()

# 画PCA分组
dat.pca <- PCA(t(exp_merge_after), graph = F)
# 提取PCA坐标数据
pca_data <- dat.pca$ind$coord[, c(1:2)]
pca_data <- as.data.frame(pca_data)
pca_data$sample <- rownames(pca_data)
group <- data.frame(sample = pca_data$sample, group = group_list_merge)
pca_data <- merge(pca_data, group, by = "sample")
colnames(pca_data)[2:3] <- c("Dim1","Dim2")
pca1 <- paste0('Dim1', ' (', round(dat.pca$eig[1, 2], 1), '%)')
pca2 <- paste0('Dim2', ' (', round(dat.pca$eig[2, 2], 1), '%)')
# 
p1 <- ggplot(pca_data, aes(x = Dim1, y = Dim2, color = group))+ 
  geom_point(aes(fill = group), size = 3, alpha = 0.7, shape = 21, position = "jitter")+
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF"))+
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF"))+
  geom_vline(xintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  labs(x = pca1, y = pca2, title = 'PCA Score Plot')+
  guides(color = "none", fill = guide_legend(title = "Group")) +
  ggprism::theme_prism()+
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12, color = 'black'),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 14),
        # legend.background = element_rect(fill = alpha("white",0)),
        legend.background = ggfun::element_roundrect(color = "black", linetype = 1, size = 0.8),
        plot.title = element_text(face = "bold", hjust = 0.50, size = 15),
        # panel.border = element_rect(color = 'black', linewidth = 0.8),
        axis.ticks = element_line(color = 'black', linewidth = 0.8))
p1
ggsave("./PCA_diff.pdf", p1, width = 8, height = 6, dpi = 600)

#### 做一个差异看看
library(limma)
design <- model.matrix( ~ 0 + group_list_merge)
colnames(design) <- levels(factor(group_list_merge, levels = c("normal", "tumor")))
rownames(design) <- colnames(exp_merge_after)
#Cancer写左，Normal写右
contrast.matrix <- makeContrasts(tumor - normal, levels = design)
contrast.matrix
#非线性最小二乘法
fit <- lmFit(exp_merge_after, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
#用经验贝叶斯调整t-test中方差的部分
fit2 <- eBayes(fit2)
#voom
nrDEG_limma_voom <- topTable(fit2, coef = 1, n = Inf)
nrDEG_limma_voom <- na.omit(nrDEG_limma_voom)
nrDEG_limma_voom_sig <- nrDEG_limma_voom %>% 
  data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(logFC != "NA") %>%
  mutate(type = case_when(logFC >= 1 & adj.P.Val < 0.05 ~ "sigUp",
                          logFC <= -1 & adj.P.Val < 0.05 ~ "sigDown",
                          .default = "nonSig"))
# 查看CHGA的logFC
nrDEG_limma_voom_sig[nrDEG_limma_voom_sig$gene_id %in% "CHGA", ]
write.csv(nrDEG_limma_voom_sig, "./nrDEG_limma_voom_sig.csv", quote = F, row.names = F)

### 
# log化
exp_merge_log <- apply(exp_merge_after[ , 1:ncol(exp_merge_after)], 2, function(x){round(log2(x+1),3)}) %>% 
  as.data.frame()
### 箱式图
library(ggplot2)
library(ggpubr)
library(ggtext)
# 提取CHGA
CHGA <- exp_merge_log %>% 
  filter(rownames(exp_merge_log) == "CHGA") %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  rename("value" = "CHGA") %>% 
  cbind(group = group$group) %>% 
  mutate(group = factor(group, levels = c("normal", "tumor")))
table(CHGA$group)
# 
p <- ggplot(CHGA, aes(x = group, y = value))+
  gghalves::geom_half_violin(aes(fill = group, color = group), side = "l", alpha = 0.8, 
                             position = position_nudge(x = -0.1, y = 0))+
  gghalves::geom_half_point(aes(fill = group, color = group), alpha = 0.5, side = "l", size = 3,
                            position = position_nudge(x = 0.1, y = 0))+
  gghalves::geom_half_boxplot(aes(fill = group, color = group), side = "r", alpha = 0.8, width = 0.25, size = 1) +
  scale_fill_manual(values = c("#1F77B4FF","#D62728FF")) +
  scale_color_manual(values = c("#1F77B4FF","#D62728FF")) +
  geom_signif(comparisons = list(c("normal", "tumor")), test = "t.test",
              tip_length = 0, size = 1, textsize = 7,
              vjust = 0.5, map_signif_level = T)+ # 用*表示显著性，*---0.05，**---0.01，***---0.001
  scale_y_continuous(labels = function(label) sprintf("%4.1f", label))+ # 设置y轴数值保留小数位数
  ggprism::theme_prism(base_size = 12) +
  labs(x = 'Group', y = 'CHGA expression-log<sub>2</sub>(count+1)', title = 'Expression from microarray profiles')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown())
p
ggsave("./half_violin_plot.pdf", p, units = "cm", width = 12, height = 13, dpi = 600)
#### 
