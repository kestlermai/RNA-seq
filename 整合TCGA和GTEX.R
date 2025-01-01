library(tidyverse)
library(data.table)
library(edgeR)
library(sva)

# 从https://xenabrowser.net/datapages/下载的数据
setwd("C:/Users/maihuanzhuo/Desktop/RNA-seq/bulk RNA-seq")
### 读取TCGA
TCGA_COAD <- fread("./TCGA/TCGA-COAD.star_counts.tsv.gz", data.table = F)
TCGA_annotation <- fread("./TCGA/gencode.v36.annotation.gtf.gene.probemap", data.table = F)
#### 读取phenotype
TCGA_pheno <- fread("./TCGA/TCGA-COAD.clinical.tsv.gz") %>%
  filter(sample %in% colnames(TCGA_COAD)) %>%
  mutate(type = ifelse(endsWith(sample, "01A"), "tumor", "normal"),
         .after = sample) %>%
  arrange(type)

## 更改TCGA表达矩阵基因名
TCGA_annotation <- TCGA_annotation %>% 
  select(id,gene)
colnames(TCGA_COAD)[1] <- "id"
TCGA_COAD <- inner_join(TCGA_annotation, TCGA_COAD) %>% 
  select(-id) %>% 
  as.data.frame()
TCGA_COAD_duplicate <- avereps(TCGA_COAD[, -1], ID = TCGA_COAD$gene) %>% 
  as.data.frame()

# order count matrix
TCGA_COAD_duplicate <- TCGA_COAD_duplicate[, TCGA_pheno$sample]

# 读取GTEX
GTEX_infor <- fread("./GTEX/probeMap_gencode.v23.annotation.gene.probemap", data.table = F)
GTEX <- fread("./GTEX/gtex_gene_expected_count.gz", data.table = F)
GTEX_phenotype <- fread("./GTEX/GTEX_phenotype.gz", data.table = F)
table(GTEX_phenotype$`_primary_site`)
## 筛选GTEx中的结直肠样本
GTEX_Colon <- GTEX %>% 
  select(sample, which(colnames(GTEX) %in% GTEX_phenotype[GTEX_phenotype$`_primary_site` == "Colon", ]$Sample))
## 更改GTEx表达矩阵基因名
GTEX_infor <- GTEX_infor %>% 
  select(id, gene)
colnames(GTEX_Colon)[1] <- "id"
## 合并数据，并删除 gene id保留gene symbol
GTEX_Colon <- inner_join(GTEX_infor, GTEX_Colon) %>% 
  select(-id) %>% 
  as.data.frame()
# 这里选择edgeR包中的函数进行取均值去重
GTEX_Colon_duplicate <- avereps(GTEX_Colon[, -1], ID = GTEX_Colon$gene) %>% 
  as.data.frame()

# devtools::install_github("wuchx101876/BioRookie")
# 可以用博主开发的R包将count值转换成TPM值，但差异分析推荐用count不建议用TPM、FPKM
# 转成TPM要把还原回去原始的count值，因为下载的数据本来就经过log2转换

## 数据合并
GTEX_Colon_duplicate$gene <- rownames(GTEX_Colon_duplicate)
TCGA_COAD_duplicate$gene <- rownames(TCGA_COAD_duplicate)
# 
combine_data <- inner_join(GTEX_Colon_duplicate, TCGA_COAD_duplicate)
combine_data_raw <- combine_data %>% 
  select(-gene) %>% 
  as.data.frame()
rownames(combine_data_raw) <- combine_data$gene

# 过滤低表达基因
combine_data_clean <- combine_data_raw[edgeR::filterByExpr(combine_data_raw, min.count = 1),]
# 取整
combine_data_integer <- apply(combine_data_clean, c(1,2), as.integer)# 包括行、列都进行取整处理

############## 整理metadata
### 处理批次效应
GTEX_Colon_phenotype <- GTEX_phenotype %>% 
  filter(`_primary_site` == "Colon", Sample %in% colnames(GTEX))
# 这里
metadata <- data.frame(sample = c(colnames(GTEX_Colon_duplicate)[-length(colnames(GTEX_Colon_duplicate))], TCGA_pheno$sample),
                       group = c(rep("normal", ncol(GTEX_Colon_duplicate) - 1), TCGA_pheno$type),
                       sex = c(GTEX_Colon_phenotype$`_gender`, TCGA_pheno$gender.demographic), 
                       batch = c(rep("GTEx", ncol(GTEX_Colon_duplicate) - 1), rep("TCGA", nrow(TCGA_pheno)))) %>%
  mutate(var = paste(batch, group, sep = "-"))
## 去批次前
# get cpm 获取标准化的矩阵
cpm_df <- cpm(y = combine_data_integer, log = T, normalized.lib.sizes = T)
# 画个PCA看看
library(FactoMineR)
library(factoextra)
dat.pca <- PCA(t(cpm_df), graph = F)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = metadata$var, 
                         palette = c("#D62728FF","#1F77B4FF","#FF7F0EFF"),
                         addEllipses = T, # Concentration ellipses
                         legend.title = "Group",
                         title = "PCA Score Plot")+
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white",0)),
        plot.title = element_text(hjust = 0.50, size = 15),
        panel.border = element_rect(color = 'black', linewidth = 1.0),
        axis.ticks = element_line(color = 'black', linewidth = 0.8))
pca_plot

# GTEX很明显有离群值，提取个体坐标和样本名称
ind_coords <- as.data.frame(dat.pca$ind$coord)
ind_coords$SampleID <- rownames(ind_coords)
outlier_sample <- ind_coords %>% 
  filter(abs(Dim.1) > 100) %>% 
  select(SampleID)
# 从cpm_df中移除离群样本
combine_data_filtered <- combine_data_integer[, !colnames(combine_data_integer) %in% outlier_sample$SampleID]
# 从metadata中移除离群样本（如果有元数据）
metadata_filtered <- metadata[!metadata$sample %in% outlier_sample$SampleID, ]
table(metadata_filtered$sex)
##### 这里性别有两个缺失，删除掉影响不大。
sex_NA <- metadata_filtered %>% 
  filter(sex == "")
metadata_filtered <- metadata_filtered %>% 
  filter(sex != "")
combine_data_filtered <- combine_data_filtered[, !colnames(combine_data_filtered) %in% sex_NA$sample]
# 重新标化
cpm_df_filtered <- cpm(y = combine_data_filtered, log = T, normalized.lib.sizes = T)
# 画图
dat.pca <- PCA(t(cpm_df_filtered), graph = F)
# 提取PCA坐标数据
pca_data <- dat.pca$ind$coord[, c(1:2)]
pca_data <- as.data.frame(pca_data)
pca_data$sample <- rownames(pca_data)
pca_data <- merge(pca_data, metadata_filtered, by = "sample")
colnames(pca_data)[2:3] <- c("Dim1","Dim2")
pca1 <- paste0('Dim1', ' (', round(dat.pca$eig[1, 2], 1), '%)')
pca2 <- paste0('Dim2', ' (', round(dat.pca$eig[2, 2], 1), '%)')
# 
p1 <- ggplot(pca_data, aes(x = Dim1, y = Dim2, color = var))+ 
  geom_point(aes(fill = var), size = 3, alpha = 0.7, shape = 21, position = "jitter")+
  stat_ellipse(aes(fill = var), geom = 'polygon', level = 0.95, alpha = 0.1, linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF"))+
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF"))+
  geom_vline(xintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  labs(x = pca1, y = pca2, title = 'Before Batch Correction')+
  guides(color = "none", fill = guide_legend(title = "Group")) +
  # theme_bw()+
  # ggpubr::theme_classic2() +
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
ggsave("./results/PCA_before_combat.pdf", p1, width = 8, height = 6, dpi = 600)

################# combat去批次
df_combat <- sva::ComBat_seq(counts = as.matrix(combine_data_filtered),
                             batch = metadata_filtered$batch,
                             group = metadata_filtered$group)
# 标化
df_combat_cpm <- cpm(y = df_combat, log = T, normalized.lib.sizes = T)
dat.pca_after <- PCA(t(df_combat_cpm), graph = F)
# 提取PCA坐标数据
pca_data_after <- dat.pca_after$ind$coord[, c(1:2)]
pca_data_after <- as.data.frame(pca_data_after)
pca_data_after$sample <- rownames(pca_data_after)
pca_data_after <- merge(pca_data_after, metadata_filtered, by = "sample")
colnames(pca_data_after)[2:3] <- c("Dim1","Dim2")
pca1 <- paste0('Dim1', ' (', round(dat.pca_after$eig[1, 2], 1), '%)')
pca2 <- paste0('Dim2', ' (', round(dat.pca_after$eig[2, 2], 1), '%)')
# 
p2 <- ggplot(pca_data_after, aes(x = Dim1, y = Dim2, color = var))+ 
  geom_point(aes(fill = var), size = 3, alpha = 0.7, shape = 21, position = "jitter")+
  stat_ellipse(aes(fill = var), geom = 'polygon', level = 0.95, alpha = 0.1, linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF"))+ # "#D62728FF","#1F77B4FF","#FF7F0EFF"
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF"))+
  geom_vline(xintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = 'black', linetype = "dashed", linewidth = 0.8) +
  labs(x = pca1, y = pca2, title = 'After Batch Correction')+
  guides(color = "none", fill = guide_legend(title = "Group")) +
  # theme_bw()+
  # ggpubr::theme_classic2() +
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
p2
ggsave("./results/PCA_after_combat.pdf", p2, width = 8, height = 6, dpi = 600)
# 合并图
p_merge <- ggpubr::ggarrange(p1, p2, common.legend = T, legend = "bottom")
p_merge
ggsave("./results/PCA_batch_merge.pdf", p_merge, width = 12, height = 6, dpi = 600)
### 保存一下
save(df_combat, metadata_filtered, file = './TCGA-GTEx.Rdata')
rm(list = ls())
# 差异分析
library(DESeq2)
load('./TCGA-GTEx.Rdata')
## 
coldata <- data.frame(row.names = colnames(df_combat),
                      condition = factor(metadata_filtered$group, levels = c("normal", "tumor")),
                      sex = factor(metadata_filtered$sex, levels = c("male", "female")))
# DESeq2输入不需要提前标化，用原始整数count matrix
df_combat <- floor(2^df_combat - 1)# xena都是以及转换成log2(x+1)
dds <- DESeqDataSetFromMatrix(countData = df_combat,
                              colData = coldata, design = ~ sex + condition)
dds <- DESeq(dds, parallel = T)
# 
DESeq2_res <- results(dds,contrast = c("condition","tumor","normal")) %>% # 指定比较的对比组condition：tumor对比normal
  data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(log2FoldChange != "NA") %>%
  mutate(type = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "sigUp",
                          log2FoldChange <= -1 & padj < 0.05 ~ "sigDown",
                          .default = "nonSig"))
table(DESeq2_res$type)
# 查看CHGA的logFC
DESeq2_res[DESeq2_res$gene_id %in% "CHGA", ]
write.csv(DESeq2_res, "./results/DESeq2_diff.csv", quote = F, row.names = F)

### 换edgeR看看，同样edgeR也是输入原始整数count matrix
DGElist <- DGEList(counts = df_combat, group = coldata$condition)
# # 过滤低质量数据
# keep <- filterByExpr(DGElist)
# DGElist <- DGElist[keep,,keep.lib.sizes = FALSE]
# TMM校正
DGElist <- calcNormFactors(DGElist)
design <- model.matrix(~ coldata$sex + coldata$condition)
# 数据离散估计 Estimating the dispersion
# estimateDisp(y, design) ：这是一个方便的包装函数，一次调用即可依次执行上述三个步骤。
# 它首先估计共同离散度，然后是趋势离散度，最后是基因特异的离散度。
DGElist <- estimateDisp(DGElist, design)
## 估计共同离散度，假设所有基因的离散度相同
# DGElist <- estimateGLMCommonDisp(DGElist, design)
## 估计趋势离散度，允许离散度随着基因表达水平的变化而变化，捕捉表达水平与离散度之间的关系
# DGElist <- estimateGLMTrendedDisp(DGElist, design)
## 估计基因特异的离散度（tagwise dispersion） ，为每个基因估计各自的离散度，允许不同基因有不同的变异性。
# DGElist <- estimateGLMTagwiseDisp(DGElist, design)
# 选上面其中之一
plotBCV(DGElist) # 数据离散分布
# 差异表达矩阵
# 拟合广义线性模型
fit <- glmFit(DGElist, design)
# 执行似然比检验
lrt <- glmLRT(fit, coef = 3)# design矩阵的第二列（即肿瘤和健康）对应的系数进行检验
topTags(lrt, n = 10)
# 返回所有结果
diff.table <- topTags(lrt, n = nrow(DGElist))$table
edgeR_diff <- diff.table %>% 
  data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(logFC != "NA") %>%
  mutate(type = case_when(logFC >= 1 & FDR < 0.05 ~ "sigUp",
                          logFC <= -1 & FDR < 0.05 ~ "sigDown",
                          .default = "nonSig"))
# 查看CHGA的logFC
edgeR_diff[edgeR_diff$gene_id %in% "CHGA", ]
write.csv(edgeR_diff, "./results/edgeR_diff.csv", quote = F, row.names = F)

######### 
# 在比较不同样本或基因之间的表达水平时，通常需要进行标准化，例如使用TPM、FPKM、RPKM或CPM等方法，
# 以消除测序深度和基因长度等因素的影响

# RPKM（Reads Per Kilobase of transcript per Million mapped reads）
# RPKM的原理是将每个基因的测量读数（reads count）标准化，以便比较不同基因的表达水平，并考虑到基因长度和总读数。
# 这种标准化方法旨在纠正测量深度和基因长度之间的差异，以便更准确地估计基因的相对表达水平。

# FPKM（Fragments Per Kilobase of transcript per Million mapped fragments）
# FPKM的原理是将每个基因的测量片段数（fragment count）标准化，以便比较不同基因的表达水平，并考虑到基因长度和总片段数。
# 和前面的RPKM一样，它也是为了纠正测量深度和基因长度之间的差异，以便更准确地估计基因的相对表达水平。
# FPKM其实是对RPKM的改进，考虑的是测量的片段数，因为在某些情况下，一个读数可能包含多个片段。

# TPM（Transcripts Per Million）
# TPM的原理是将每个基因或转录本的测量次数（counts）标准化，以便比较不同基因或转录本的表达水平，
# 同时考虑到总的测量次数。TPM的核心思想是将每个基因或转录本的表达水平标准化为每百万个映射转录本的表达水平，
# 从而减小测量深度和样本之间的影响。

# CPM（Counts Per Million）
# 除了前面的 RPKM、FPKM、TPM等，CPM 也是一种较为常见的衡量基因或转录本表达水平的标准化方法。
# CPM的原理是将每个基因或转录本的测量次数（counts）标准化，以便比较不同基因或转录本的表达水平，并考虑到总的测量次数。
# 它将每个基因或转录本的表达水平标准化为每百万次测量次数的表达水平，从而减小测量深度和样本之间的影响。

# RPKM：适合基因的相对表达水平比较；适用于单端测序数据。
# FPKM：适合基因的相对表达水平比较；适用于双端测序数据。
# TPM：基因的相对表达水平比较；适用于不同样本之间的比较。
# CPM：适合基因的相对表达水平比较；适用于简单分析任务。

### 这里考虑用TPM
# 把count转换成TPM需要gene length
# https://mp.weixin.qq.com/s?__biz=MzU4NjU4ODQ2MQ==&mid=2247490627&idx=1&sn=f23242af5baa6cd6c07ff3558d65c97b&chksm=fdf85401ca8fdd17017107413706c4c259c91385ddc17e0778539d327f2f1d49f0f5072ea688&scene=21#wechat_redirect
# 嫌麻烦直接用别人封装好的函数了，应该是推荐是https://mp.weixin.qq.com/s/XXpc0IBa2SJlFp5iTH1OtA
# 可以直接用IOBR包直接跑出来
# devtools::install_github("IOBR/IOBR")
library(IOBR)
df_tpm <- count2tpm(countMat = df_combat, idType = "SYMBOL", org = "hsa")
# log化
df_tpm_log <- apply(df_tpm[ , 1:ncol(df_tpm)], 2, function(x){round(log2(x+1),3)}) %>% 
  as.data.frame()

### 箱式图
library(ggplot2)
library(ggpubr)
library(ggtext)
# 提取CHGA
CHGA <- df_tpm_log %>% 
  filter(rownames(df_tpm_log) == "CHGA") %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::rename("value" = "CHGA") %>%
  cbind(group = metadata_filtered$group) %>% 
  mutate(group = factor(group, levels = c("normal", "tumor")))

table(CHGA$group)
# normal  tumor 
# 362    453 
p <- ggplot(CHGA, aes(x = group, y = value))+
  gghalves::geom_half_violin(aes(fill = group, color = group), side = "l", alpha = 0.8, 
                             position = position_nudge(x = -0.1, y = 0))+
  gghalves::geom_half_point(aes(fill = group, color = group), alpha = 0.5, side = "l", size = 1.5,
                            position = position_nudge(x = 0.1, y = 0))+
  gghalves::geom_half_boxplot(aes(fill = group, color = group), side = "r", alpha = 0.8, width = 0.2, size = 1) +
  scale_fill_manual(values = c("#1F77B4FF","#D62728FF")) +
  scale_color_manual(values = c("#1F77B4FF","#D62728FF")) +
  geom_signif(comparisons = list(c("normal", "tumor")), test = "wilcox.test",
              tip_length = 0, size = 1, textsize = 7,
              vjust = 0.5, map_signif_level = T)+ # 用*表示显著性，*---0.05，**---0.01，***---0.001
  scale_y_continuous(labels = function(label) sprintf("%4.1f", label))+ # 设置y轴数值保留小数位数
  ggprism::theme_prism(base_size = 10) +
  labs(x = 'Group', y = 'CHGA Expression-log<sub>2</sub>(TPM+1)', title = 'Expression from TCGA and GTEx mRNA profiles')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown())
p
ggsave("./results/half_violin_plot.pdf", p, units = "cm", width = 12, height = 13, dpi = 600)

############# 生存分析
### 读取TCGA
TCGA_COAD <- fread("./TCGA/TCGA-COAD.star_counts.tsv.gz", data.table = F)
TCGA_annotation <- fread("./TCGA/gencode.v36.annotation.gtf.gene.probemap", data.table = F)
#### 读取phenotype
TCGA_pheno <- fread("./TCGA/TCGA-COAD.clinical.tsv.gz") %>%
  filter(sample %in% colnames(TCGA_COAD)) %>%
  mutate(type = ifelse(endsWith(sample, "01A"), "tumor", "normal"),
         .after = sample) %>%
  arrange(type)
colnames(TCGA_pheno)
## 更改TCGA表达矩阵基因名
TCGA_annotation <- TCGA_annotation %>% 
  select(id,gene)
colnames(TCGA_COAD)[1] <- "id"
TCGA_COAD <- inner_join(TCGA_annotation, TCGA_COAD) %>% 
  select(-id) %>% 
  as.data.frame()
TCGA_COAD_duplicate <- avereps(TCGA_COAD[, -1], ID = TCGA_COAD$gene) %>% 
  as.data.frame()
# 提取生存资料
col_group <- c("sample",#样本ID
               "gender.demographic",
               "days_to_death.demographic",#死亡时间
               "days_to_last_follow_up.diagnoses",#最后随访时间
               "vital_status.demographic",#生or死
               "ajcc_pathologic_stage.diagnoses",#临床分期
               "age_at_diagnosis.diagnoses"#初始病例诊断年龄
)
### 只筛选tumor数据
tumor_pheno <- TCGA_pheno %>% 
  filter(type == "tumor") %>% 
  select(col_group) %>% 
  # 把生存状态转换为0(Alive)和1(Dead):
  mutate(event = ifelse(vital_status.demographic == "Alive", 0, 1),
         # 死亡时间和最后随访时间缺失值替换为0
         days_to_death.demographic = ifelse(is.na(days_to_death.demographic), 0, days_to_death.demographic),
         days_to_last_follow_up.diagnoses = ifelse(is.na(days_to_last_follow_up.diagnoses), 0, days_to_last_follow_up.diagnoses),
         # 计算生存时间（死亡时间+最后随访时间）---探讨生存情况，而不是诊断年龄+死亡时间---这个是计算生存率
         OS_time = as.numeric(days_to_death.demographic) + as.numeric(days_to_last_follow_up.diagnoses)) %>% 
  # 因为病人个人原因，而导致的随访时间过短，这些数据没有意义要去掉；
  filter(OS_time >= 30)# 去掉生存时间小于30天的样本

### 筛选肿瘤样本表达矩阵
survival_df <- TCGA_COAD_duplicate[, tumor_pheno$sample]
survival_df <- floor(2^survival_df - 1)# xena都是以及转换成log2(x+1)
## 在bulk RNA-seq中只有差异分析用到raw count，其他用TPM、FPKM、CPM这些标化数据
survival_df_tpm <- count2tpm(countMat = survival_df, idType = "SYMBOL", org = "hsa")
# log化
survival_df_tpm_log <- apply(survival_df_tpm[ , 1:ncol(survival_df_tpm)], 2, function(x){round(log2(x+1),3)}) %>% 
  as.data.frame()
### 用中位数/均数划分表达量高低：常见用中位数
#as.integer转化为整数
tumor_pheno$group <- ifelse(as.integer(survival_df_tpm_log["CHGA",]) > 
                              median(as.integer(survival_df_tpm_log["CHGA",])), "high", "low")
## ok 开始生存分析
library(survminer)
library(survival)
#用survival包survfit函数拟合生存曲线
fit <- survfit(Surv(OS_time, event) ~ group, data = tumor_pheno)
fit
#cox回归
cox_fit <- coxph(Surv(OS_time, event) ~ group, data = tumor_pheno)
cox_fit
# 
pdf("./results/survival_curve.pdf", width = 7, height = 7)
ggsurvplot(fit, data = tumor_pheno,
           censor.shape = "|", censor.size = 4,
           conf.int = TRUE,
           conf.int.style = "ribbon",
           conf.int.alpha = 0.2,
           pval = "logrank p-value = 0.11 \n HR (low) = 0.48 \n p-value (HR) = 0.11",
           palette = "lancet",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           legend = "top",
           legend.labs = c("High","Low"),
           xlab = "OS_time (days)",
           ylab = "Survival probablity",
           title = "Survival curves",
           break.x.by = 1000,
           break.y.by = 0.2,
           risk.table = TRUE,
           risk.table.col = "strata",
           risk.table.height = 0.2,
           risk.table.y.text = T)
dev.off()

##### 画一个不同癌症分期的小提琴图
tumor_pheno$value <- t(survival_df_tpm_log["CHGA",])
tumor_pheno_stage <- tumor_pheno %>% 
  filter(ajcc_pathologic_stage.diagnoses != "") %>% 
  # 去除I、II、III或IV后面的ABC，这样感觉太麻烦了而且写的不对，改成对存在大写的ABC都去除
  # mutate(stage = str_replace(ajcc_pathologic_stage.diagnoses, "(Stage [IV|I{1,3}])A|B|C", "\\1"),
  #        stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))
  mutate(stage = str_replace(ajcc_pathologic_stage.diagnoses, "[ABC]", ""), # 去掉大写的 A, B, C
         stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))
table(tumor_pheno_stage$stage)
# combn函数用于生成给定向量中所有可能的组合
comp <- combn(c("Stage I", "Stage II", "Stage III", "Stage IV"), 2)
my_comparisons <- list()
for(i in 1:ncol(comp)){
  my_comparisons[[i]] <- comp[,i]
}
my_comparisons
# 
p <- ggpubr::ggviolin(tumor_pheno_stage, x = "stage", y = "value",
                      fill = "stage", add = "boxplot", add.params = list(fill = "white"))+
  # stat_compare_means(comparisons = my_comparisons) +
  stat_pwc(method = "t.test", tip.length = 0.03, label.size = 5, size = 0.8,
           bracket.nudge.y = 0.2, label = "p.signif") +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  ggprism::theme_prism(base_size = 14) +
  labs(x = 'Group', y = 'Expression-log<sub>2</sub>(TPM+1)', title = 'CHGA from TCGA mRNA profile')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown())
p
ggsave("./results/stage_violin_plot.pdf", p, units = "cm", width = 16, height = 14, dpi = 600)
# 
p <- ggplot(tumor_pheno_stage, aes(x = stage, y = value))+
  gghalves::geom_half_violin(aes(fill = stage, color = stage), side = "l", alpha = 0.8, 
                             position = position_nudge(x = -0.05, y = 0))+
  gghalves::geom_half_point(aes(fill = stage, color = stage), alpha = 0.6, side = "l", size = 2,
                            position = position_nudge(x = 0.05, y = 0))+
  gghalves::geom_half_boxplot(aes(fill = stage, color = stage), side = "r", alpha = 0.8, width = 0.2, size = 1) +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF"))+
  stat_pwc(method = "t.test", tip.length = 0.03, label.size = 5, size = 0.8,
           bracket.nudge.y = 0.1, label = "p.signif") +
  scale_y_continuous(labels = function(label) sprintf("%4.1f", label))+ # 设置y轴数值保留小数位数
  ggprism::theme_prism(base_size = 14) +
  labs(x = 'Group', y = 'Expression-log<sub>2</sub>(TPM+1)', title = 'CHGA from TCGA mRNA profile')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown())
p
ggsave("./results/stage_half_violin_plot.pdf", p, units = "cm", width = 16, height = 14, dpi = 600)
# 

##### 从https://ibdmdb.org/下载的转录组数据
df_raw <- fread("./ibdmdb/host_tx_counts.tsv.gz", data.table = F) # 没有log化直接用
metadata <- read.csv("./ibdmdb/hmp2_metadata_2018-08-20.csv")
# 现在处理一下metadata
metadata_filtered <- metadata %>% 
  filter(data_type == "host_transcriptomics") %>% 
  mutate(group = case_when(diagnosis == "UC" ~ "IBD",
                           diagnosis == "CD" ~ "IBD",
                           diagnosis == "nonIBD" ~ "Control")) %>% 
  select(External.ID, group) %>% 
  dplyr::rename(sample = External.ID)
### 处理转录组表达矩阵
mat <- df_raw %>% 
  column_to_rownames(var = "Transcript")
### reads起码是百万级别
column_totals <- data.frame(Sample = colnames(mat), Total = colSums(mat[, sapply(mat, is.numeric)]))
filtered_samples <- column_totals[column_totals$Total > 1000000, ]
mat_filter <- mat[,filtered_samples$Sample]

### 转成TPM
library(IOBR)
mat_tpm <- IOBR::count2tpm(countMat = mat_filter, idType = "SYMBOL", org = "hsa")
mat_tpm_log <- apply(mat_tpm[ , 1:ncol(mat_tpm)], 2, function(x){round(log2(x+1),3)}) %>% 
  as.data.frame()

### 比较组间差异
library(ggplot2)
library(ggpubr)
library(ggtext)
# 提取CHGA
CHGA <- mat_tpm_log %>% 
  filter(rownames(mat_tpm_log) == "CHGA") %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::rename("value" = "CHGA") %>% 
  inner_join(metadata_filtered, by = "sample") %>% 
  mutate(group = factor(group, levels = c("Control", "IBD")))
table(CHGA$group)
# Control     IBD 
# 50     196
# 箱式图
p <- ggplot(CHGA, aes(x = group, y = value))+
  gghalves::geom_half_violin(aes(fill = group, color = group), side = "l", alpha = 0.8, 
                             position = position_nudge(x = -0.1, y = 0))+
  gghalves::geom_half_point(aes(fill = group, color = group), alpha = 0.5, side = "l", size = 1.5,
                            position = position_nudge(x = 0.1, y = 0))+
  gghalves::geom_half_boxplot(aes(fill = group, color = group), side = "r", alpha = 0.8, width = 0.2, size = 1) +
  scale_fill_manual(values = c("#1F77B4FF","#D62728FF")) +
  scale_color_manual(values = c("#1F77B4FF","#D62728FF")) +
  geom_signif(comparisons = list(c("Control", "IBD")), test = "wilcox.test",
              tip_length = 0, size = 1, textsize = 7,
              vjust = 0.5, map_signif_level = T)+ # 用*表示显著性，*---0.05，**---0.01，***---0.001
  scale_y_continuous(labels = function(label) sprintf("%4.1f", label))+ # 设置y轴数值保留小数位数
  ggprism::theme_prism(base_size = 12) +
  labs(x = 'Group', y = 'CHGA Expression-log<sub>2</sub>(TPM+1)', title = 'RNA profile from IBDMDB project')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown(face = "bold"))
p
ggsave("./ibdmdb/half_violin_plot.pdf", p, units = "cm", width = 12, height = 13, dpi = 600)

#### 整个矩阵包含转录本即包括mRNA、rRNA、tRNA、ncRNA（非编码）
# 注释文件用的是GRCh38.p14，https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz
### gencode.v47.annotation.gtf.gz感觉太新
# 这里考虑还是用GRCh38.p13的gencode.v36.annotation.gtf.gz
library(rtracklayer)
anno <- import("./ibdmdb/gencode.v36.annotation.gtf.gz")
anno <- as.data.frame(anno)#将文件转换为数据框格式
pc_anno <- anno %>% 
  filter(gene_type == "protein_coding", type == "gene")
protein_coding_genes <- pc_anno$gene_name
mat_anno <- mat_filter[rownames(mat_filter) %in% protein_coding_genes, ]
# 去除行名中的小数点及其后的部分，这是基因的版本
rownames(mat_anno) <- sub("\\..*$", "", rownames(mat_anno))
# 如果有重复可以取平均

# 差异分析
library(DESeq2)
coldata <- data.frame(row.names = CHGA$sample,
                      condition = factor(CHGA$group, levels = c("Control", "IBD")))
# DESeq2输入不需要提前标化，用原始整数count matrix，要求整数，而且不是log转化过的（注意xena的数据）
dds <- DESeqDataSetFromMatrix(countData = mat_anno[, CHGA$sample], colData = coldata, design = ~ condition)
dds <- DESeq(dds, parallel = T)
# estimating size factors
# Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
# every gene contains at least one zero, cannot compute log geometric means
# 只有在所有样品中表达都不为0的基因才能用来计算
# dds <- DESeq(dds, parallel = T, sfType = "iterate") # iterate，这样算的太慢了，从头开始做一下质控吧，现在就能跑了
### 
# 
DESeq2_res <- results(dds, contrast = c("condition","IBD","Control")) %>% # 指定比较的对比组condition：tumor对比normal
  data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(log2FoldChange != "NA") %>%
  mutate(type = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "sigUp",
                          log2FoldChange <= -1 & padj < 0.05 ~ "sigDown",
                          .default = "nonSig"))
table(DESeq2_res$type)
# 查看CHGA的logFC
DESeq2_res[DESeq2_res$gene_id %in% "CHGA", ]
write.csv(DESeq2_res, "./ibdmdb/DESeq2_diff.csv", quote = F, row.names = F)

#### 单基因GSEA
library(org.Sc.sgd.db)
library(clusterProfiler)
# 
group <- ifelse(as.numeric(mat_anno["CHGA", CHGA$sample]) > median(as.numeric(mat_anno["CHGA", CHGA$sample])), 
                "High", "Low")
group <- factor(group, levels = c("Low", "High"))
### 重新跑一次DESeq2
coldata <- data.frame(row.names = CHGA$sample, condition = group)
# DESeq2输入不需要提前标化，用原始整数count matrix，要求整数，而且不是log转化过的（注意xena的数据）
dds <- DESeqDataSetFromMatrix(countData = mat_anno[, CHGA$sample], colData = coldata, design = ~ condition)
dds <- DESeq(dds, parallel = T)
# 
DESeq2_res <- results(dds, contrast = c("condition","High","Low")) %>% # 指定比较的对比组condition：tumor对比normal
  data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  filter(log2FoldChange != "NA") %>%
  mutate(type = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "sigUp",
                          log2FoldChange <= -1 & padj < 0.05 ~ "sigDown",
                          .default = "nonSig"))
table(DESeq2_res$type)
# 查看CHGA的logFC
DESeq2_res[DESeq2_res$gene_id %in% "CHGA", ]
write.csv(DESeq2_res, "./ibdmdb/DESeq2_diff-GESA.csv", quote = F, row.names = F)
#
genelist <- DESeq2_res %>% 
  dplyr::select(gene_id, log2FoldChange) %>% 
  arrange(desc(log2FoldChange))

### 太久没跑富集分析了，都忘了富集分析之前需要把ENTREZID
ENTREZID <- bitr(geneID = genelist$gene_id, 
                 fromType = "SYMBOL", # 需要转换ID类型
                 toType = "ENTREZID", # 转换成的ID类型
                 OrgDb = "org.Hs.eg.db") # 对应的物种，小鼠的是org.Mm.eg.db
# 使用merge合并
gesa_list <- merge(genelist, ENTREZID, by.x = "gene_id", by.y = "SYMBOL", all = F)
# 按照logFC降序排序
gesa_list_sort <- gesa_list[order(gesa_list$log2FoldChange, decreasing = T),]
# 把foldchange按照从大到小提取出来
gene_list <- gesa_list_sort$log2FoldChange
# 给上面提取的foldchange加上对应上ENTREZID
names(gene_list) <- gesa_list_sort$ENTREZID
head(gene_list)
save(gene_list, file = "./ibdmdb/GSEA-gene_list.Rdata")
load("./ibdmdb/GSEA-gene_list.Rdata")
#开始富集分析
kegg_ges <- gseKEGG(geneList = gene_list,
                    organism = 'hsa',
                    nPerm = 1000,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    eps = 0)
kegg_result <- as.data.frame(kegg_ges)
write.csv(kegg_result, "./ibdmdb/GSEA-KEGG.csv",  quote = F, row.names = F)
# devtools::install_github("junjunlab/GseaVis", force = T)
library(GseaVis)
#GO enrich
BP_ges <- gseGO(geneList = gene_list,
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 0)
BP_ges_result <- BP_ges@result
#### 
gseaNb(object = kegg_ges, 
       geneSetID = kegg_ges@result$ID[1],
       subPlot = 3,
       curveCol = "#996699",
       lineSize = 1.2, 
       geneSize = 5,
       pCol = 'black', 
       addPval = T,
       pvalX = 0.95,
       pvalY = 0.8,
       pDigit = 3,
       nesDigit = 3)
ggsave(filename = "./ibdmdb/GSEA-KEGG-IL17.pdf", width = 18, height = 14, units = "cm", dpi = 600)

### 突然想到如果是UC和CD都能看见下降是不是更好说明呢
df_raw <- fread("./ibdmdb/host_tx_counts.tsv.gz", data.table = F) # 没有log化直接用
metadata <- read.csv("./ibdmdb/hmp2_metadata_2018-08-20.csv")
# 现在处理一下metadata
metadata_filtered <- metadata %>% 
  filter(data_type == "host_transcriptomics") %>% 
  mutate(group = case_when(diagnosis == "UC" ~ "UC",
                           diagnosis == "CD" ~ "CD",
                           diagnosis == "nonIBD" ~ "Control")) %>% 
  dplyr::select(External.ID, group) %>% 
  dplyr::rename(sample = External.ID)
### 处理转录组表达矩阵
mat <- df_raw %>% 
  column_to_rownames(var = "Transcript")
### reads起码是百万级别
column_totals <- data.frame(Sample = colnames(mat), Total = colSums(mat[, sapply(mat, is.numeric)]))
filtered_samples <- column_totals[column_totals$Total > 1000000, ]
mat_filter <- mat[,filtered_samples$Sample]
### 转成TPM
library(IOBR)
mat_tpm <- IOBR::count2tpm(countMat = mat_filter, idType = "SYMBOL", org = "hsa")
mat_tpm_log <- apply(mat_tpm[ , 1:ncol(mat_tpm)], 2, function(x){round(log2(x+1),3)}) %>% 
  as.data.frame()

### 比较组间差异
library(ggplot2)
library(ggpubr)
library(ggtext)
# 提取CHGA
CHGA <- mat_tpm_log %>% 
  filter(rownames(mat_tpm_log) == "CHGA") %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  dplyr::rename("value" = "CHGA") %>% 
  inner_join(metadata_filtered, by = "sample") %>% 
  mutate(group = factor(group, levels = c("Control", "CD", "UC")))
table(CHGA$group)
# Control      CD      UC 
# 50     123      73 

p <- ggplot(CHGA, aes(x = group, y = value))+
  gghalves::geom_half_violin(aes(fill = group, color = group), side = "l", alpha = 0.8, 
                             position = position_nudge(x = -0.1, y = 0))+
  gghalves::geom_half_point(aes(fill = group, color = group), alpha = 0.5, side = "l", size = 1.5,
                            position = position_nudge(x = 0.1, y = 0))+
  gghalves::geom_half_boxplot(aes(fill = group, color = group), side = "r", alpha = 0.8, width = 0.2, size = 1) +
  scale_fill_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF")) +
  scale_color_manual(values = c("#D62728FF","#1F77B4FF","#FF7F0EFF")) +
  geom_signif(comparisons = list(c("Control", "CD"),
                                 c("Control", "UC"),
                                 c("CD", "UC")), test = "t.test",
              tip_length = 0.02, size = 1, textsize = 6,
              vjust = 0, y_position = c(11, 12, 11), extend_line = -0.01,
              map_signif_level = T)+ # 用*表示显著性，*---0.05，**---0.01，***---0.001
  scale_y_continuous(labels = function(label) sprintf("%4.1f", label))+ # 设置y轴数值保留小数位数
  ggprism::theme_prism(base_size = 11) +
  labs(x = 'Group', y = 'CHGA Expression-log<sub>2</sub>(TPM+1)', title = 'Transcriptional Profile from the IBDMDB Project')+
  guides(color = 'none', fill = 'none') +
  theme(axis.title.y = element_markdown(face = "bold"))
p
ggsave("./ibdmdb/half_violin_plot-CD-UC.pdf", p, units = "cm", width = 13, height = 14, dpi = 600)

#### 将ENTREZID转成SYMBOL
library(org.Sc.sgd.db)
library(clusterProfiler)
list <- c(3627, 6356, 4318, 3458, 6347, 2920, 7124, 2921, 3565, 3383, 24, 112744, 1437, 6280, 6372, 
          2919, 6374, 5743, 6354, 8061, 27189, 6279, 6278, 3553, 3596, 1673, 3605, 4586, 4312, 4314, 
          4322, 3569, 1440)
SYMBOL.ID <- bitr(geneID = list, 
                  fromType = "ENTREZID", # 需要转换ID类型
                  toType = "SYMBOL", # 转换成的ID类型
                  OrgDb = "org.Hs.eg.db") # 对应的物种，小鼠的是org.Mm.eg.db
write.csv(SYMBOL.ID, "./ibdmdb/IL17_SYMBOL_ID.csv",  quote = F, row.names = F)
