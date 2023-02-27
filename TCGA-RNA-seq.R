setwd("C:/Users/11146/Desktop/R包/TCGA")
#TCGA表达矩阵TCGA表达矩阵中的ENSG_ID包含了编码基因、非编码基因以及假基因
#对mRNA或lncRNA单独研究，这时候就必须将原始表达矩阵做拆分
#制作包含geneID和gene type信息的mRNA和lncRNA清单
#下载GDC官方注释ENSG_ID文件 https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

#读取GDC官方的V22版本ENS注释文件：
gcodev22 <- read.table("gencode.gene.info.v22.tsv",header = T,row.names = 1,check.names = F)

#行名（ENSG_ID）、gene_name列以及gene_type列即为我们所需
#gene_name用于ID转换，gene_type用于提取mRNA和lncRNA，现在我们提取这三列，并合并为一个数据框。

#提取gene name列：
gname <- subset(gcodev22,select = c(1,6))#提取第1,6列
gname <- cbind(gene_stable_ID=row.names(gname), gname)#把gname的行名设置为第一列，列名为gene_stable_ID
row.names(gname)=NULL#去掉原行名
head(gname)

#拆分为仅包含mRNA和lncRNA的两个独立数据框。
#mRNA在gene_type中对应的是protein_coding，但lncRNA在gene_type中不止对应lncRNA，相关信息我们可以从gencode数据库获取。
#根据gene_type分别提取lncRNA和mRNA对应的ENSID和geneID：
#lncRNA的type来自https://www.gencodegenes.org/pages/biotypes.html的lncRNA部分
lncRNA <- c(
  '3prime_overlapping_ncRNA',
  'antisense',
  'bidirectional_promoter_lncRNA',
  'lincRNA',
  'macro_lncRNA',
  'non_coding',
  'processed_transcript',
  'sense_intronic',
  'sense_overlapping')

mRNA <- c("protein_coding")

#分别生成mRNA、lncRNA的ENSID、geneID和type对应list：
lncRNA_list <- gname[gname$gene_type %in% lncRNA,]#%in%逻辑运算符，检查左边元素是否出现在右边中，若有，返T，数据框保存返回T的数据
mRNA_list <- gname[gname$gene_type %in% mRNA,]
head(lncRNA_list)
head(mRNA_list)

#保存TCGA原始矩阵
save(lncRNA_list,mRNA_list,file = c("gene_ID_list.Rdata"))

load("gene_ID_list.Rdata")

#读取Xena下载的KICH(嫌色细胞癌)表达量矩阵(counts)：
#或者用TCGAbiolink下载，后续不用log转换
exp <- read.table("TCGA-KICH.htseq_counts.tsv",header = T,row.names = 1,check.names = F)
exp <- as.matrix(exp)#将数据转换成矩阵
exp[1:6,1:6]#查看前6个

#将表达量矩阵还原为原始counts
exp <- 2^exp-1
#取整数
exp <- round(exp)
exp[1:6,1:6]
dim(exp)#获取行=基因数，列=样本数

#去掉最后5行：
exp <- exp[-(length(exp[,1]):(length(exp[,1])-4)),]

#数据过滤，TCGA为了保证geneID在不同癌症的一致性，导致很多样本的表达量为0
exp_filtered <- exp[apply(exp, 1, function(x) sum(x > 1) > 89*0.5), ]#保留75%的样本基因数大于1；过滤的标准可以是保留50%，exp大于10的基因
dim(exp_filtered)

#分别取表达矩阵和mRNA/lncRNA_list的ENSID交集：
#mRNA交集：
mRNA <- intersect(rownames(exp_filtered),mRNA_list$gene_stable_ID)

#lncRNA交集：
lncRNA <- intersect(rownames(exp_filtered),lncRNA_list$gene_stable_ID)

#查看过滤后的mRNA和lncRNA数量(交集部分)：
length(mRNA)
length(lncRNA)

#mRNA表达矩阵：从表达矩阵中提取交集部分(mRNA)：
mRNA_exp <- exp_filtered[mRNA,]

#lncRNA表达矩阵：从表达矩阵中提取交集部分(lncRNA)：
lncRNA_exp <- exp_filtered[lncRNA,]

dim(mRNA_exp)
dim(lncRNA_exp)#行数和交集gene数相等，确认无误，表达矩阵拆分完成
#保存仅包含mRNA和lncRNA的两个独立表达矩阵
save(mRNA_exp,lncRNA_exp,mRNA_list,lncRNA_list,file = c('mRNA_lncRNA_ENSID_exp.Rdata'))

#样本分组，对TCGA样本名称的第14、15位判断样本属于“tumor” or “normal”
#TCGA-02-0001-01C-01D-0182-01；0-9为“tumor”，10-19为“normal”
#提取样本（列名）中第14-15位的字符串并查看
library(stringr)
table(str_sub(colnames(exp_filtered),14,15))

#不是所有数据都有tumor和normal的样本
#创建分组（normal和tumor），并将其转化为因子指定顺序：
group_list = ifelse(as.numeric(str_sub(colnames(exp_filtered),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#保存矩阵跟分组
save(exp_filtered,group_list,file = c("TCGA_KICH_exp.Rdata"))

#主流用raw_count进行差异分析，RPKM、FPKM已经过时
#DESeq2要求用raw_count，然后要做的是多个样本间同一特征比较前的均一化，并不适合样本内比较，不同基因的表达水平
#RPKM、FPKM、TPM是样本内的所有特征做均一化（基因长度均一化、测序长度均一化），并不适合样本间差异分析
#基因长度归一化的原因：同样表达水平下，某基因长度越长，对应得到的reads数越多。归一化后「同一样本的不同基因表达水平」间具有可比性；
#测序深度归一化的原因：不同样本的mRNA reads总量有高有底，归一化后「不同样本的基因表达水平」间具有可比性。
#DESeq2差异分析
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
#导入整理好的表达矩阵
load("TCGA_KICH_exp.Rdata")

#创建表达量矩阵样本名和分组一一对应的数据框（colData）
colData <- data.frame(row.names =colnames(exp_filtered),condition=group_list)
head(colData)

#用函数DESeqDataSetFromMatrix()构建DESeqDataSet对象（dds），包括存储输入值、计算中间值、标准化处理等
dds <- DESeqDataSetFromMatrix(countData = exp_filtered,colData = colData, design = ~ condition)

#构建广义线性模型
dds <- DESeq(dds)

#因子设置tumor在前，normal在后
res <- results(dds,contrast = c("condition",rev(levels(group_list))))
res <- as.data.frame(res)
head(res)
#保存差异分析数据：
save(dds,res,file = c("TCGA_KICH_DESeq2.Rdata"))

#DESeq2包中三种数据标准化方法：
#blind=T则每个样本之间自动计算标准化因子
vsd <- vst(dds, blind=FALSE)#variance stabilizing transformations
rld <- rlog(dds, blind=FALSE)#regularized logarithm-rlog，常用但计算时间长
ntd <- normTransform(dds)#常规log2(n + 1)

rld <- as.data.frame(assay(rld))#获取标化后的矩阵
#保存矩阵
save(rld,file=c("TCGA_KICH_count_transformation.Rdata"))

load("TCGA_KICH_count_transformation.Rdata")

#简单可视化
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("tinyarray")
library(tinyarray)
library(dplyr)
library(pheatmap)

#PCA,用标化的rld
KICH_PCA <- draw_pca(rld,group_list)
KICH_PCA

#volano plot，火山图用count
#生成显著上下调gene标签列：
res$group <- case_when(res$log2FoldChange > 1 & res$pvalue < 0.05 ~ "Up",
                       res$log2FoldChange < -1 & res$pvalue < 0.05 ~ "Down",
                       abs(res$log2FoldChange) <= 1 ~ "None",
                       res$pvalue >= 0.05 ~ "None")
head(res)

#heatmap，热图也一样用count
#选取标准化后平均表达量Top10绘制：
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]

KICH_heatmap <- draw_heatmap(rld[select,],
                            group_list,
                            scale_before = F,
                            legend = T,
                            annotation_legend = T,
                            color_an = c("pink", "#8DA0CB"))
KICH_heatmap

#选取所有显著上下调gene绘制：
select2 <- rownames(res)[res$group !="None"]

KICH_heatmap2 <- draw_heatmap(rld[select2,],
                              group_list,
                              scale_before = F,
                              legend = T,
                              annotation_legend = T,
                              color_an = c("#97cb00", "#ff938a"),
                              color = (grDevices::colorRampPalette(c("orange", "white","purple")))(100))
KICH_heatmap2

#保存图片data：
save(KICH_PCA,KICH_volcano,KICH_heatmap,KICH_heatmap2,file = c("TCGA_KICH_DESeq2_plot.Rdata"))

#ensembl gene ID对应的gene name
#geneID转换
load("mRNA_lncRNA_ENSID_exp.Rdata")
#行名顺序匹配
#按照mRNA矩阵行名的顺序来匹配mRNA_list中gene_stable_ID的顺序，对mRNA_list进行顺序重排并取交集：
a <- rownames(mRNA_exp)
b <- mRNA_list$gene_stable_ID
mRNA_ID_match <- mRNA_list[match(a,b),]

head(mRNA_ID_match)
dim(mRNA_ID_match)

#lncRNA同理
c <- rownames(lncRNA_exp)
d <- lncRNA_list$gene_stable_ID
lncRNA_ID_match <- lncRNA_list[match(c,d),]

head(lncRNA_ID_match)
dim(lncRNA_ID_match)

#行名不能重复，去除重复geneID
#去掉mRNA表达矩阵中有着相同gene name的不同ENSG_ID：
#duplicated函数用来判断
dp_mRNA <- duplicated(mRNA_ID_match$gene_name)
table(dp_mRNA)

#但我们需要的是不重复的gene（逻辑值为FALSE），因此做一个逻辑非运算（即给出相反逻辑值，将不重复的gene返回为TRUE）
n = !dp_mRNA
#根据逻辑值取子集，筛选mRNA表达矩阵中返回逻辑值为TRUE的行以及全部列
dp_mRNA_exp <- mRNA_exp[n,]
#在完成顺序匹配的mRNA_ID_match中进行相同规则筛选
dp_mRNA_ID_match <- mRNA_ID_match[n,]
dim(dp_mRNA_exp)
dim(dp_mRNA_ID_match)#检查去重

#修改表达矩阵中行名的ENSID为gene name(已完成顺序匹配和去重复)
rownames(dp_mRNA_exp) <- dp_mRNA_ID_match$gene_name
dp_mRNA_exp[1:6,1:4]#gene ID转换完成

dp_lncRNA <- duplicated(lncRNA_ID_match$gene_name)
table(dp_lncRNA)
n=!dp_lncRNA
dp_lncRNA_exp <- lncRNA_exp[n,]
dp_lncRNA_ID_match <- lncRNA_ID_match[n,]
dim(dp_lncRNA_exp)
dim(dp_lncRNA_ID_match)

rownames(dp_lncRNA_exp) <- dp_lncRNA_ID_match$gene_name
dp_lncRNA_exp[1:6,1:4]

#保存
save(dp_mRNA_exp,dp_lncRNA_exp,file = c("mRNA_lncRNA_geneID_exp.Rdata"))

#读取临床信息表格和载入表达矩阵：
cl_KICH <- data.table::fread("TCGA-KICH.GDC_phenotype.tsv")
#拆分完的mRNA和lncRNA表达矩阵并已经转换geneID
load("mRNA_lncRNA_geneID_exp.Rdata")

cl_KICH <- as.data.frame(cl_KICH)
#用mRNA表达矩阵
mRNA_exp_clean <- dp_mRNA_exp
dim(cl_KICH)

#提取临床信息表格中所需要的列名(需要根据自己的癌症背景知识挑选)；
##常规选择：如病人ID、生死状态、性别、人种、临床分期等
##创建所需列名的分组文件：
col_group <- c(
  "submitter_id.samples",#样本ID
  "submitter_id",#病人ID
  "days_to_death.demographic",#死亡时间
  "days_to_last_follow_up.diagnoses",#最后随访时间
  "vital_status.demographic",#生or死
  "tumor_stage.diagnoses",#临床分期
  "age_at_initial_pathologic_diagnosis"#初始病例诊断年龄
  #"more......性别、人种等等"
)

#筛选列名：
cl_KICH2 <- cl_KICH[,col_group]
dim(cl_KICH)#临床信息184个样本
dim(mRNA_exp_clean)#表达矩阵89个样本

#首先，我们需要去掉表达矩阵中的正常样本（样本编号第14-15位字符串≥10为normal）：
##查看tumor和normal样本数量：
table(str_sub(colnames(mRNA_exp_clean),14,15))#表达矩阵中tumor样本65个，normal样本24个

#创建tumor和normal的分组文件：
sample_group <- ifelse(as.numeric(str_sub(colnames(mRNA_exp_clean),14,15)) < 10,"tumor","normal")

#筛选肿瘤样本：
mRNA_exp_clean2 <- mRNA_exp_clean[,sample_group == c("tumor")]
dim(mRNA_exp_clean2)#只剩65个肿瘤样本了

#将临床表格中样本行名的顺序按照表达矩阵列名进行匹配：
cl_KICH3 <- cl_KICH2[match(colnames(mRNA_exp_clean2),cl_KICH2$submitter_id.samples),]
dim(cl_KICH3)#匹配完也只剩65个肿瘤样本

#检查顺序是否一致以及一一对应：
head(colnames(mRNA_exp_clean2))
head(cl_KICH3$submitter_id.samples)
identical(colnames(mRNA_exp_clean2),cl_KICH3$submitter_id.samples)


#把临床信息表格第一列作为行名：
rownames(cl_KICH3) = cl_KICH3[,1]
cl_KICH3 <- cl_KICH3[,-1]
View(cl_KICH3)
View(mRNA_exp_clean2)

######把生存状态转换为0(Alive)和1(Dead):
event <- cl_KICH3$vital_status.demographic
event <- ifelse(event =="Alive",0,1)
#添加event列
cl_KICH3$event <- event
table(cl_KICH3$event)

#先把所有缺失的NA替换为0：
colnames(cl_KICH3)
#把死亡时间和最后随访时间缺失值替换为0
for (i in 1:ncol(cl_KICH3[,2&3])){
  cl_KICH3[,i][is.na(cl_KICH3[,i])] <- 0
}

#把记录时间有问题的样本整行剔除：
cl_KICH4 <- cl_KICH3[!rownames(cl_KICH3) %in% c("TCGA-KN-8430-01A","TCGA-KL-8343-01A"),]
#同样删除表达矩阵的单列
mRNA_exp_clean3 <- mRNA_exp_clean2[,!colnames(mRNA_exp_clean2) %in% c("TCGA-KN-8430-01A","TCGA-KL-8343-01A")]
#检查是否一一匹配：
identical(colnames(mRNA_exp_clean3),rownames(cl_KICH4))
#根据死亡时间和最后随访时间计算OS_time
cl_KICH4$OS_time <- (as.numeric(cl_KICH4[,2])+as.numeric(cl_KICH4[,3]))#两者相加得出生存时间

#诊断KICH年龄分组
age <- cl_KICH4$age_at_initial_pathologic_diagnosis
#查看病理诊断年龄区间：
range(age)

#把年龄分为青年(17-34)、中年(35-54)、老年(55-86)三组：
age <- ifelse(age<=34,"young",
              ifelse(age<=54,"middle","old"))
cl_KICH4$age <- age
table(cl_KICH4$age)

######按某基因的高低表达分组：
#选某个关键基因
gene <- c("HSF1")
#用中位数/均数划分表达量高低：常见用中位数
#as.integer转化为整数
cl_KICH4$gene <- ifelse(as.integer(mRNA_exp_clean3[gene,]) > median(as.integer(mRNA_exp_clean3[gene,])), "high", "low")
cl_KICH4$gene <- ifelse(as.integer(mRNA_exp_clean3[gene,]) > mean(as.integer(mRNA_exp_clean3[gene,])), "high", "low")
cl_KICH4[1:6,7:10]

######保存整理完成临床信息和表达矩阵
save(cl_KICH4,mRNA_exp_clean3,file = c("survival_data.Rdata"))

library(survminer)
library(survival)
load("survival_data.Rdata")
colnames(cl_KICH4)
#探究stage分期对KICH生存率的影响
#用survival包survfit函数拟合生存曲线
fit_stage <- survfit(Surv(OS_time, event) ~ tumor_stage.diagnoses, data = cl_KICH4)
#绘图
ggsurvplot(fit_stage, cl_KICH4,
           cencor.shape = "|", cencor.size = 4,#删失点形状，default“+”
           conf.int = T, conf.int.style = "ribbon", #置信区间类型，默认"ribbon",可选"step(虚线)"
           conf.int.alpha = 0.2,#置信区间不透明度调节
           pval = T,#p值
           palette = "lancet",
           ggtheme =theme_bw(),
           legend = "right",
           legend.labs = c("stage Ⅰ","stage Ⅱ","stage Ⅲ","stage Ⅳ"),#
           xlab = "OS_time(days)",
           ylab = "Survival Probablity",
           title = "Survival Curves",
           break.x.by = 1000,
           break.y.by = 0.2,
           #add.all = T,添加总生存曲线，即所有病人不分组
           surv.median.line = "hv")#添加中位生存时间线，“hv”、“h”、“v”，v为绘制垂直线，h为绘制水平线

#添加风险表或删失事件图：
p <- ggsurvplot(fit_stage, cl_KICH4,
           cencor.shape = "|", cencor.size = 4,#删失点形状，default“+”
           conf.int = T, conf.int.style = "ribbon", #置信区间类型，默认"ribbon",可选"step(虚线)"
           conf.int.alpha = 0.2,#置信区间不透明度调节
           pval = T,#p值
           palette = "lancet",
           ggtheme =theme_bw(),
           legend = "right",
           legend.labs = c("stage Ⅰ","stage Ⅱ","stage Ⅲ","stage Ⅳ"),#
           xlab = "OS_time(days)",
           ylab = "Survival Probablity",
           title = "Survival Curves",
           break.x.by = 1000,
           break.y.by = 0.2,
           #add.all = T,添加总生存曲线，即所有病人不分组
           surv.median.line = "hv",
           risk.table = TRUE,#风险表添加
           risk.table.col = "strata",#风险表颜色跟随
           risk.table.height = 0.2,#生存表高度占据画幅百分比(区间0-1，1为只显示生存表)
           risk.table.y.text = FALSE, #隐藏风险表y轴标签
           ncensor.plot = TRUE, #删失事件图绘制
           ncensor.plot.height = 0.15 #删失事件图高度占据画幅百分比(区间0-1，同上)
           )#添加中位生存时间线，“hv”、“h”、“v”，v为绘制垂直线，h为绘制水平线

# add method to grid.draw
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("survial curve.png", p, width = 8, height = 10, dpi = 600)

#探究age对KICH生存率的影响
fit_age <- survfit(Surv(OS_time, event) ~ age, data = cl_KICH4)

ggsurvplot(
  fit_age,
  data = cl_KICH4,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  palette = "lancet",
  #surv.median.line = "hv",
  ggtheme =theme_bw(),
  legend = "right",
  legend.labs = c("Middle","Old","Young"),
  xlab = "OS_time(days)",
  ylab = "Survival Probablity",
  title = "Survival Curves",
  break.x.by = 1000,
  break.y.by = 0.2,
  risk.table = TRUE,#风险表添加
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)

#探究单个gene对KICH生存率的影响
#上面选的关键基因是HSF1
fit_gene_HSF1 <- surv_fit(Surv(OS_time,event)~gene,cl_KICH4)

ggsurvplot(
  fit_gene_HSF1,
  data = cl_KICH4,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  palette = "lancet",
  #surv.median.line = "hv",#某组的中位生存时间还没达到，所以无法添加中位生存时间线
  ggtheme =theme_bw(),
  legend = "top",
  legend.labs = c("High","Low"),
  xlab = "OS_time(days)",
  ylab = "Survival probablity",
  title = "Survival curves",
  break.x.by = 1000,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)

#HSF1对生存率无显著影响，换成ARHGAP33
rownames(mRNA_exp_clean3)
gene <- c("ARHGAP33")
cl_KICH4$gene <- ifelse(mRNA_exp_clean3[gene,] > median(mRNA_exp_clean3[gene,]),"High","Low")
cl_KICH4$gene[1:10]

fit_gene_ARHGAP33 <- surv_fit(Surv(OS_time,event)~gene,cl_KICH4)

ggsurvplot(
  fit_gene_ARHGAP33,
  data = cl_KICH4,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  palette = "lancet",
  #surv.median.line = "hv",#某组的中位生存时间还没达到，所以无法添加中位生存时间线
  ggtheme =theme_bw(),
  legend = "top",
  legend.labs = c("High","Low"),
  xlab = "OS_time(days)",
  ylab = "Survival probablity",
  title = "Survival curves",
  break.x.by = 1000,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)
#ARHGAP33高表达会影响KICH患者生存率

#cox回归
cox_age <- coxph(Surv(OS_time, event) ~ age_at_initial_pathologic_diagnosis, data = cl_KICH4)
cox_age
#看看ARHGAP33在cox回归有没有意义
gene <- c("ARHGAP33")
cl_KICH4$gene <- ifelse(mRNA_exp_clean3[gene,] > median(mRNA_exp_clean3[gene,]),"High","Low")

cox_gene_ARHGAP33 <- coxph(Surv(OS_time, event) ~ gene, data = cl_KICH4)
summary(cox_gene_ARHGAP33)#没意义，除了p值，exp（coef）即HR风险比也要关注

#使用exp连续变量
cl_KICH4$gene <- mRNA_exp_clean3[gene,]
cox_gene_ARHGAP33 <- coxph(Surv(OS_time, event) ~ gene, data = cl_KICH4)
summary(cox_gene_ARHGAP33)

ggsurvplot(survfit(cox_gene_ARHGAP33), cl_KICH4)#用survfit计算cox的生存函数

#批量对所有基因进行生存分析,批量计算HR值等
#apply函数(x矩形，1=对行操作，2=对列操作，3=对每个元素进行操作)
cox_results <- apply(
  mRNA_exp_clean3,1,function(x){
    cl_KICH4$gene <- ifelse(x>median(x),"High","Low")
    cox_gene <- coxph(Surv(OS_time, event) ~ gene, data = cl_KICH4)
    beta <- coef(cox_gene) #回归系数，即coef值
    se <- sqrt(diag(vcov(cox_gene))) #标准误
    HR <- exp(beta) #HR(Hazard Ratio/风险比)值计算
    HRse <- HR * se
    tmp <- round(cbind(coef = beta,
                            se = se,
                            z = beta/se, #z值，即Wald test统计量
                            p = 1 - pchisq((beta/se)^2, 1),#显著性P值
                            HR = HR, #风险比
                            HRse = HRse,
                            HRz = (HR - 1) / HRse,
                            HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                            #95%CI：
                            ower_95 = exp(beta - qnorm(.975, 0, 1) * se),
                            upper_95 = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    return(tmp["geneLow",])
  }
)
cox_results <- t(cox_results)
head(cox_results)[,1:3]

#查看统计学上有显著差异的基因数量：
table(cox_results[,4] < 0.05)

#提取有显著性差异的基因列表：
cox_results_sig5 <- cox_results[cox_results[,4]<0.05,]
head(cox_results_sig5)[,1:3]


#批量计算log rank p值
#用survdiff函数计算多条生存曲线之间差异，得到卡方chisq值
logrank_p <- apply(
  mRNA_exp_clean3,1,function(x){ 
    cl_KICH4$gene <- ifelse(x>median(x),"High","Low")
    diff <- survdiff(Surv(OS_time, event) ~ gene, data=cl_KICH4)#默认logrank检验
    p <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    #pchisq 函数是用来计算卡方分布的累积分布函数值 (cumulative distribution function)
    #length(diff$n) - 1 是卡方分布的自由度，表示可以自由变化的变量数
    return(p)
  }
)
#按p值从小到大排序,并将其转化为数据框：
logrank_p_order <- as.data.frame(sort(logrank_p),header = T)

#提取显著影响的基因--对应行名
logrank_p_sig5 <- rownames(logrank_p_order)[logrank_p_order<0.05]
logrank_p_sig1 <- rownames(logrank_p_order)[logrank_p_order<0.01]
head(logrank_p_sig1)
length(logrank_p_sig1)

##采取批量绘图方式，将有显著性差异Top6 基因绘制生存曲线图：
#lapply函数是对list每个向量进行操作，apply是对行或列进行操作
genes_top6 <- c(logrank_p_sig1[1:6])
p_logrank_top6 <- lapply(genes_top6, function(x){
  cl_KICH4$gene <- ifelse(mRNA_exp_clean3[x,] > median(mRNA_exp_clean3[x,]),"High","Low")
  fit <- survfit(Surv(OS_time, event) ~ gene, data=cl_KICH4)
  ggsurvplot(
    fit,
    data = cl_KICH4,
    censor.shape="|", censor.size = 4,
    conf.int = TRUE,
    conf.int.style = "ribbon",
    conf.int.alpha = 0.2,
    pval = TRUE,
    palette = "lancet",#npg
    surv.median.line = "hv",
    ggtheme =theme_bw(),
    legend = "top",
    legend.labs = c("High","Low"),
    xlab = "OS_time(days)",
    ylab = "Survival Probablity",
    title = paste0(x," Signature"),
    break.x.by = 1000,
    break.y.by = 0.2,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.2,
    risk.table.y.text = FALSE
  )
}
)
#拼图：
p_logrank_top6 <- arrange_ggsurvplots(
  p_logrank_top6,
  print = T,
  title = "Survival Curves",
  ncol = 3, 
  nrow = 2,
  risk.table.height = 0.25
  )
ggsave("logrank_top6.pdf", p_logrank_top6, width = 15, height = 10, dpi = 600)

#还可以和logrank检验结果的p值取交集：
#提取显著差异的基因名
cox_results_sig_gene <- rownames(cox_results_sig)
head(cox_results_sig_gene)

#intersect函数提取两个数组共有的元素
int_p5 <- intersect(logrank_p_sig5,cox_results_sig_gene)
int_p1 <- intersect(logrank_p_sig1,cox_results_sig_gene)
head(int_p5)
length(int_p5)
head(int_p1)
length(int_p1)

colnames(cl_KICH4)

surv_int_gene <- lapply(int_p1,function(x){
  cl_KICH4$gene <- ifelse(mRNA_exp_clean3[x,] > median(mRNA_exp_clean3[x,]),"High","Low")
  fit_int_gene <- survfit(Surv(OS_time, event) ~ gene, data = cl_KICH4)
  ggsurvplot(
    fit_int_gene,
    data = cl_KICH4,
    censor.shape="|", censor.size = 4,
    conf.int = TRUE,
    conf.int.style = "ribbon",
    conf.int.alpha = 0.2,
    pval = TRUE,
    palette = "lancet",#npg
    surv.median.line = "hv",
    ggtheme =theme_bw(),
    legend = "top",
    legend.labs = c("High","Low"),
    xlab = "OS_time(days)",
    ylab = "Survival Probablity",
    title = paste0(x," Signature"),
    break.x.by = 1000,
    break.y.by = 0.2,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.1,
    risk.table.y.text = FALSE
  )
}
)
surv_int_gene <- arrange_ggsurvplots(
  surv_int_gene,
  print = T,
  title = "Survival Curves",
  ncol = 7, 
  nrow = 3,
  risk.table.height = 0.25
)
ggsave("surv_int_gene.pdf", surv_int_gene, width = 30, height = 20, dpi = 600)

save(int_p1,file = c("int_gene.Rdata"))

####ensemble_ID转换成gene_symbol
load("TCGA_KICH_DESeq2.Rdata")
head(res)
#去除ensemble_ID的版本号
library(stringr)
res$ensembl_gene_id=unlist(str_split(row.names(res),"[.]",simplify=T))[,1]
res[1:3,7]

#去除后三行没有ensemble_ID的数据
res <- res[1:(nrow(res)-3),]#前面去除了后面就不用去除

#对基因进行注释-获取gene_symbol，用bioMart对ensembl_id转换成gene_symbol
#BiocManager::install("biomaRt", force = T)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) 
save(mart,file = "mart.Rdata")
load("mart.Rdata")

hg_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), 
                   filters= 'ensembl_gene_id', 
                   values = res$ensembl_gene_id, 
                   mart = mart)

#xy是两个表合并的列名，两个表数据匹配合并，并保留两个表的所有列
res <- merge(res, hg_symbols, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

#检查是否有重复基因
any(duplicated(colnames(res)))

#将空值转成缺失值，然后去除没有注释的基因，去除4271个基因
res[res==""] <- NA
res <- na.omit(res)

save(res,file = c("res_symbol.Rdata"))

#KEGG富集分析
library(dplyr)
library(org.Hs.eg.db)#人类基因注释包
library(clusterProfiler)#富集分析包括id转换、go注释
library(ggplot2)
library(RColorBrewer)
#添加上下调基因分组标签：
res$group <- case_when(res$log2FoldChange > 2 & res$pvalue < 0.05 ~ "Up",
                       res$log2FoldChange < -2 & res$pvalue < 0.05 ~ "Down",
                       abs(res$log2FoldChange) <= 2 ~ "None",
                       res$pvalue >= 0.05 ~ "None")
head(res)

#分别筛选上调基因、下调基因或所有差异基因（上调+下调）：
up <- res$hgnc_symbol[res$group=="Up"]#差异上调
down <- res$hgnc_symbol[res$group=="Down"]#差异下调
diff <- c(up,down)#所有差异基因
head(up)
head(down)

#用clusterProfiler包进行gene_symbol转换ENTREZID，或者直接ensembl_id转成ENTREZID
##使用函数bitr(基于org.Hs.eg.db包)：
columns(org.Hs.eg.db)
#up：
up_entrez <- bitr(up,
                  fromType = "SYMBOL",#现有的ID类型 "SYMBOL"或者"ENSEMBL"
                  toType = "ENTREZID",#需转换的ID类型
                  OrgDb = "org.Hs.eg.db")
#2.13% fail to map
#down（同上）：
down_entrez <- bitr(down,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
#0.26% fail to map
#diff（同上）：
diff_entrez <- bitr(diff,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
#1.04% fail to map
head(diff_entrez)
save(diff_entrez,file = c("diff_entrez.Rdata"))
#小部分没有mapping到很正常

#KEGG富集分析，我们以总差异基因(diff_entrez)为例：
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa",#物种，Homo sapiens (human)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)


KEGG_result <- KEGG_diff@result
#保存富集结果：
save(KEGG_diff,KEGG_result,file = c("KEGG_diff.Rdata"))

#探索/调出选定的KEGG通路
#红色为富集到该通路的差异基因
browseKEGG(KEGG_diff, "hsa04610")

load("KEGG_diff.Rdata")
#富集可视化
library(enrichplot)
#条形图
barplot(
  KEGG_diff,
  x = "Count", #or "GeneRatio"
  color = "pvalue", #or "p.adjust" and "qvalue"
  showCategory = 20,#显示前top20
  font.size = 12,
  title = "KEGG enrichment barplot",
  label_format = 30 #超过30个字符串换行
)

#气泡图
dotplot(
  KEGG_diff,
  x = "GeneRatio",
  color = "p.adjust",
  title = "Top 20 of Pathway Enrichment",
  showCategory = 20,
  label_format = 30
)

#用-log10p或者-log10q来绘图
#先提取富集结果表前Top20：
KEGG_top20 <- KEGG_result[1:20,]
#指定绘图顺序（转换为因子）：
KEGG_top20$pathway <- factor(KEGG_top20$Description,levels = rev(KEGG_top20$Description))

#Top20富集数目条形图：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

p <- ggplot(data = KEGG_top20,
            aes(x = Count,
                y = pathway,
                fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdPu",direction = 1) +
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "KEGG enrichment barplot") + mytheme
p

#Top20显著富集条形图：
p1 <- ggplot(data = KEGG_top20,
             aes(x = -log10(pvalue),
                 y = pathway,
                 fill = Count)) +
  scale_fill_distiller(palette = "Blues",direction = 1) +
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "-log10(pvalue)",
       y = "pathway",
       title = "KEGG enrichment barplot") +mytheme
p1

#Top20显著性气泡图：
#将pathway按照p值排列：
p2 <- ggplot(data = KEGG_top20,
             aes(x = Count,
                 y = pathway))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+ # 气泡大小及颜色设置
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Count") +
  mytheme
p2

#将pathway按照gene number（Count）数排列：
p3 <- ggplot(data = KEGG_top20,
             aes(x = Count,
                 y = reorder(pathway,Count)))+ # 用reorder将pathway按照Count重排
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Count") +
  mytheme
p3

#展现Rich Factor
#富集因子（Rich Factor）计算：

#Rich Factor = GenRatio/BgRatio
top20_rf <- apply(KEGG_top20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  KEGG_top20_rf <- round(GeneRatio/BgRatio,2)
  KEGG_top20_rf
})
head(top20_rf)
KEGG_top20$Rich_Factor <- top20_rf

#富集因子版显著性气泡图绘制：
p3 <- ggplot(data = KEGG_top20,
             aes(x = Rich_Factor, # X轴用富集因子来映射
                 y = pathway))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Rich Factor",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Gene Number") +
  mytheme
p3

##GO富集分析
#在GO富集分析中有三个Ontology
#分子功能MF(Molecular Function)、细胞组分CC(Cellular Component)及生物过程BP(Biological Process)
load("diff_entrez.Rdata")

library(dplyr)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)

GO_MF_diff <- enrichGO(gene = diff_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       ont = "MF", #GO分支，"BP"（生物学过程）、"MF"（分子功能）和"CC"（细胞组分）。或者“all”合并
                       pAdjustMethod = "BH", #多重假设检验矫正方法
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = T) #是否将gene ID映射到gene name

GO_MF_result <- GO_MF_diff@result


GO_BP_diff <- enrichGO(gene = diff_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = T) 
GO_BP_result <- GO_BP_diff@result

GO_CC_diff <- enrichGO(gene = diff_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = T) 
GO_CC_result <- GO_CC_diff@result

GO_all_diff <- enrichGO(gene = diff_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       ont = "all", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = T) 
GO_all_result <- GO_all_diff@result

##保存GO富集结果：
save(GO_MF_diff,GO_CC_diff,GO_BP_diff,GO_all_diff,file = c("GO_diff.Rdata"))
save(GO_MF_result,GO_BP_result,GO_CC_result,GO_all_result,file = c("GO_diff_result.Rdata"))

load("GO_diff_result.Rdata")
#####快速可视化探索
#BiocManager::install("topGO")
library(topGO)
library(enrichplot)
library(ggplot2)

#GO有向无环图绘制
goplot(GO_MF_diff) #来自enrichplot
plotGOgraph(GO_MF_diff) #来自topGO

#GO富集条形图：
barplot(GO_MF_diff,x = "Count", #or "GeneRatio"
  color = "pvalue", #or "p.adjust" and "qvalue"
  showCategory = 20,#显示前top20(enrichResult按照p值排序)
  font.size = 12,title = "Cellular Component enrichment barplot",
  label_format = 30) #超过30个字符串换行

#GO富集气泡图：
dotplot(GO_MF_diff,x = "GeneRatio",color = "p.adjust",
  title = "Top 20 of GO CC terms Enrichment",showCategory = 20,
  label_format = 30)

#富集网络图:
#pairwise_termsim函数是计算两者相似性，两个GO terms之间存在overlap即表明两者相关性，overlap越高相关性越高
edo <- pairwise_termsim(GO_MF_diff)
emapplot(edo,layout = "kk", #布局形式
         showCategory = 30) #展示GO terms的数量

####使用ggplot2进行可视化:
#取前top20，并简化命名：
MF_20 <- GO_MF_result[1:20,]
CC_20 <- GO_CC_result[1:20,]
BP_20 <- GO_BP_result[1:20,]

#自定义主题
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))
#在MF的Description中存在过长字符串，我们将长度超过50的部分用...代替：
MF_20$Description <- str_trunc(MF_20$Description,width = 50,side = "right")
MF_20$Description

#指定绘图顺序（转换为因子）：
MF_20$term <- factor(MF_20$Description,levels = rev(MF_20$Description))
CC_20$term <- factor(CC_20$Description,levels = rev(CC_20$Description))
BP_20$term <- factor(BP_20$Description,levels = rev(BP_20$Description))

#GO富集柱形图：
GO_bar <- function(x){
  y <- get(x)
  ggplot(data = y,aes(x = Count,y = term,fill = -log10(pvalue))) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) + #label换行，部分term描述太长
    geom_bar(stat = "identity",width = 0.8) +
    labs(x = "Gene Number",y = "Description",title = paste0(x," of GO enrichment barplot")) +
    theme_bw() +
    mytheme
}


p_bar_MF <- GO_bar("MF_20")+scale_fill_distiller(palette = "Blues",direction = 1)
p_bar_BP <- GO_bar("BP_20")+scale_fill_distiller(palette = "Reds",direction = 1)
p_bar_CC <- GO_bar("CC_20")+scale_fill_distiller(palette = "Greens",direction = 1)

#计算富集因子(Rich Factor):富集基因集中显著富集的基因数目，与全部参考基因集中对应基因数目的比值；
#富集因子的阈值大于1.5或2被认为是显著富集的
#用apply计算每一行（1）的富集基因数/总基因数=Rich Factor
#MF:
rf<- apply(MF_20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))#对GeneRatio这一列字符串转化成数字类型eval和parse函数一般同时使用
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)#保留两位小数
  RF
})
MF_20$Rich_Factor <- rf

rf<- apply(CC_20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
CC_20$Rich_Factor <- rf

rf<- apply(BP_20,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
BP_20$Rich_Factor <- rf

#绘制GO富集气泡图：
GO_dot <- function(x){
  y = get(x)
  ggplot(data = y,aes(x = Rich_Factor,y = term)) +
    geom_point(aes(size = Count,color = -log10(pvalue))) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) +
    labs(x = "Rich Factor",y = "Description",title = paste0(x,"of GO enrichment Dotplot"), 
         size = "Gene Number") + 
    theme_bw()+
    mytheme
}

#查看色板描述。
brewer.pal.info

p_dot_MF <- GO_dot("MF_20") + scale_color_distiller(palette = "YlGnBu",direction = 1)
p_dot_BP <- GO_dot("BP_20") + scale_color_distiller(palette = "YlOrRd",direction = 1)
p_dot_CC <- GO_dot("CC_20") + scale_color_distiller(palette = "YlGn",direction = 1)


#将三个ontology拉通取top30(按照p值排序)绘图：
all_result <- arrange(GO_all_result,pvalue) #默认升序

#取top30：
all_30 <- all_result[1:30,]

#指定绘图顺序（转换为因子）：
all_30$term <- factor(all_30$Description,levels = rev(all_30$Description))

#自定义y轴标签颜色(区分不同ontology)：
col_function <- function(x){
  col <- rep("black", length(x))
  BP <- which(x %in% c("BP"))
  CC <- which(x %in% c("CC"))
  MF <- which(x %in% c("MF"))
  col[BP] <- "#e6576a"
  col[CC] <- "#6D3BC6"
  col[MF] <- "#6aefdd"
  col
}
y_text_color <- col_function(all_30$ONTOLOGY)

#富集条形图
p_bar_all <- ggplot(data = all_30, aes(x = Count,y = term, fill = -log10(pvalue))) +
  geom_bar(stat = "identity",width = 0.8) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1)+
  labs(x = "Gene Number", y = "Description", title = "TOP30 GO enrichment barplot") +
  theme_bw() +
  mytheme+
  theme(axis.text.y = element_text(color = y_text_color))
p_bar_all

#计算rf
rf<- apply(all_30,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
all_30$Rich_Factor <- rf

#富集气泡图
p_dot_all <- ggplot(data = all_30, aes(x = Rich_Factor, y = term)) +
  geom_point(aes(size = Count,color = -log10(pvalue))) +
  scale_color_distiller(palette = "YlOrRd", direction = 1)+
  labs(x = "Rich Factor", y = "Description", title = "TOP30 GO enrichment Dotplot", 
       size = "Gene Number") + 
  theme_bw()+
  mytheme+
  theme(axis.text.y = element_text(color = y_text_color))
p_dot_all

#合并plot
library(patchwork)
p_bar_merge <- p_bar_MF + p_bar_BP + p_bar_CC
p_dot_merge <- p_dot_MF + p_dot_BP + p_dot_CC
p_all_merge <- p_bar_all + p_dot_all

ggsave("p_bar_merge.png", p_bar_merge, width = 20, height = 5, dpi = 600)
ggsave("p_dot_merge.png", p_dot_merge, width = 20, height = 5, dpi = 600)
ggsave("p_all_merge.png", p_all_merge, width = 15, height = 5, dpi = 600)


#GSEA富集分析
library(org.Hs.eg.db)
library(clusterProfiler)

load("res_symbol.Rdata")

##：symbol转entrez ID：
entrez <- bitr(res$hgnc_symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
#0.52% fail to map
head(entrez)

#genelist格式：entrez ID+log2fc
genelist <- res$log2FoldChange
names(genelist) <- res$hgnc_symbol
head(genelist)

#将genelist表里的symbol转换成entrez
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：
KEGG_ges <- gseKEGG(geneList = genelist,
  organism = "hsa",#"hsa"-人类,不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0)

KEGG_ges_result <- KEGG_ges@result

#GSEA可视化
##山峦图：
library(ggridges)
library(ggplot2)
library(enrichplot)

p_ridge <- ridgeplot(KEGG_ges,showCategory = 15,fill = "p.adjust",decreasing  = T)
p_ridge

p_dot <- dotplot(KEGG_ges)
p_dot


##enrichmentScore折线图：
p_score <- gseaplot(KEGG_ges,geneSetID = 1, #取第一个pathway绘图
               by = "runningScore", #or “preranked”, “all”
               title = KEGG_ges$Description[1])
p_score


##rank分布图：
p_rank <- gseaplot(KEGG_ges,geneSetID = 1,by = "preranked",title = KEGG_ges$Description[1])
p_rank


##同时显示：
p_KEGG_all <- gseaplot(KEGG_ges,
               geneSetID = 1,
               by = "all",
               title = KEGG_ges$Description[1])
p_KEGG_all

#gseaplot2函数:
##单个基因集展示：
p_KEGG_gseaplot2 <- gseaplot2(KEGG_ges,geneSetID = 1,
                color = "red",rel_heights = c(1.5, 0.5, 1),#子图高度
                subplots = 1:3,#显示哪些子图
                pvalue_table = F,#是否显示pvalue表
                title = KEGG_ges$Description[1],ES_geom = "line")#"dot"将线转换为点
p_KEGG_gseaplot2


#同时展示多个基因集：
p_KEGG_gseaplot2_merge <- gseaplot2(KEGG_ges,geneSetID = 1:3,#或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = c("#a195fb", "#ff8696", "#bdff83"),pvalue_table = T,ES_geom = "line")
p_KEGG_gseaplot2_merge


##GSEA_GO富集分析：
GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Hs.eg.db,
                ont = "CC", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges_result <- GO_ges@result

p_GO_ges_all <- gseaplot(GO_ges,
                  geneSetID = 1,
                  by = "all",
                  title = GO_ges$Description[1])
p_GO_ges_all

#在富集结果表中将entrez转换为symbol：
#setReadable函数进行ID转换
KEGG_ges_set <- setReadable(KEGG_ges,
                            OrgDb = org.Hs.eg.db,
                            keyType="ENTREZID")

KEGG_ges_set@result$core_enrichment[1]

#GSEA_Reactome富集分析:
#基于Reactome数据库中的代谢、信号传导、基因表达调控等通路信息，对基因集进行富集分析，进一步了解基因集的生物学功能和作用
#按照基因表达量从高到低排序，并在基因集中计算通路富集得分，最后使用Kolmogorov-Smirnov检验来判断每个通路是否富集
#BiocManager::install("ReactomePA")
library(ReactomePA)
ret_ges <- gsePathway(genelist,
                      organism = "human",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      eps = 0)
ret_ges_result <- ret_ges@result


#用DO(Disease Ontology)数据库进行GSEA,DOSE包
##GSEA_DO(Disease Ontology)富集分析:
library(DOSE)
DO_ges <- gseDO(genelist,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = FALSE,
                eps = 0)
DO_ges_result <- DO_ges@result

#MSigDB（GSEA提供的基因集数据库）和msigdf包配合使用
#devtools::install_github("ToledoEM/msigdf")
library(msigdf)

#提取C2注释(human)：
library(dplyr)
c2 <- msigdf.human %>%
  filter(category_code == "c2") %>% select(geneset, symbol) %>% as.data.frame

head(c2)

#重新准备genelist文件：
genelist <- res$log2FoldChange
names(genelist) <- res$hgnc_symbol
head(genelist)

#将genelist表里的symbol转换成entrez
genelist2 <- genelist[names(genelist) %in% entrez[,1]]
#将genelist按照log2FC值从高到低进行排序：
genelist2 <- sort(genelist2,decreasing = T)
head(genelist2)


c2_ges <- GSEA(genelist2,
               TERM2GENE = c2,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = FALSE,
               eps = 0)
c2_ges_result <- c2_ges@result


##关于计算富集因子（Rich Factor）和富集倍数（Fold Enrichment）的问题
#富集因子Rich Factor = 差异基因中富集到该通路的基因数（count）/背景基因中富集到该通路的基因数（BgRatio分子）
#富集倍数Fold Enrichment =GeneRatio/BgRatio

load("KEGG_diff.Rdata")
#计算Rich Factor（富集因子）：
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#计算Fold Enrichment（富集倍数）：
KEGG_diff2 <- mutate(KEGG_diff, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]

###Y叔公众号提到
#FoldEnrichment = (k/n) / (M/N)，RichFactor = k/M，
#n和N数字是不变的，所以等同于说FoldEnrichment = RichFactor * C，C = N/n 是一个常数


##富集火柴棒图
mytheme <- theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
                 legend.title = element_text(size = 13),legend.text = element_text(size = 11))

p <- ggplot(GO_CC_diff2, showCategory = 20, aes(x = RichFactor,y = reorder(Description, RichFactor))) + #将Description按照RichFactor进行排序
  geom_segment(aes(xend=0, yend = Description),size = 1.1,color = "grey") +
  geom_point(aes(color=-log10(pvalue), size = Count)) +
  scale_color_distiller(palette = "Spectral",direction = -1)+
  scale_size_continuous(range=c(2, 8)) + #气泡大小范围
  theme_minimal() +
  labs(x = "Rich factor", y = NULL, title = "GO CC terms Enrichment") + mytheme
p


#富集基因网络图
#用于展示感兴趣的通路所富集到的差异基因有哪些及上下调情况（logFC值），还能够看出哪些基因同时出现在多个通路（重叠关系）

#先制作差异基因的genelist(需要差异分析结果表格中的log2FoldChange值)：
genelist <- res$log2FoldChange
names(genelist) <- res$hgnc_symbol
length(genelist)#全部基因，16345个

##筛选差异基因：
genelist <- genelist[names(genelist) %in% diff]
length(genelist)
head(genelist)

cnetplot(KEGG_diff2,showCategory = 5,#显示top5 pathway
         categorySize = "pvalue",
         foldChange = genelist)

#自定义展示感兴趣的通路：
interest <- c("Butanoate metabolism","Histidine metabolism",
              "Proximal tubule bicarbonate reclamation","Propanoate metabolism",
              "Primary bile acid biosynthesis")

cnetplot(KEGG_diff2, showCategory = interest,
         categorySize="pvalue",
         foldChange=genelist)
#可以用layout参数调整网络图布局
#可选'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.

#另一种形式：
interest2 <- c("Butanoate metabolism","Histidine metabolism","Propanoate metabolism")

cnetplot(KEGG_diff2,
         showCategory = interest2,
         foldChange=genelist,
         circular = TRUE,
         colorEdge = TRUE)

#heatmap
#用热图的形式展示感兴趣通路所富集到的差异基因及上下调情况。
heatplot(KEGG_diff2,
         foldChange = genelist,
         showCategory = interest)

