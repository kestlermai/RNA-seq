rm(list = ls())
setwd("C:/Users/11146/Desktop/R包/免疫浸润")
options(stringsAsFactors = F)

library(tidyverse)
load("KICH_exp_tpm.Rdata")
#Cibersort混池矩阵要求:
#非负数、无缺失值、第一列为gene symbol;数据不要取log;RNA-seq测序，使用FPKM和TPM标准化数据；
#Affymetrix芯片，使用MAS5或RMA标准化数据；Illumina的Beadchip和Agilent的 single color按照limma描述处理；
#表达矩阵整理：把原本的行名（gene symbol）变为第一列；
exp <- rownames_to_column(KICH_exp_tpm)
colnames(exp)[1] <- c("Gene symbol")
write.table(exp,file = "KIRC_exp_tpm.txt",row.names = F,quote = F,sep = "\t") 

source("CIBERSORT.R")
# 首先读取免疫细胞的浸润比例LM22和RNA-seq表达矩阵
#https://doi.org/10.1038/nmeth.3337-LM22源文件
library(e1071)
library(parallel)
library(preprocessCore)
#cibersort免疫细胞以及脚本得放在同一个文件夹中
exp <- "KIRC_exp_tpm.txt"
LM22 <- "LM22.txt"

setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE33463")
load("exp.Rdata")
source("CIBERSORT.R")
exp <- cbind("gene symbol" = rownames(exp), exp)
write.table(exp,file = "exp-cibersort.txt",row.names = F,quote = F,sep = "\t") 
exp <- "exp-cibersort.txt"
LM22 <- "LM22.txt"

res <- CIBERSORT(LM22, exp, perm = 1000, #排列数量大于100计算P值
                QN = T)#数据归一化

##过滤丰度全是0的免疫细胞类型：
resm <- as.data.frame(res[,1:22])#后三列（P-value、Correlation、RMSE）去掉，绘图时暂时不需要用到
choose <- apply(resm,2,function(x) {sum(x) > 0})
table(choose)
resm <- resm[,choose]
save(res,resm,file = c("cibersort_results.Rdata"))
load("cibersort_results.Rdata")
resm <- t(resm)

library(pheatmap)
library(stringr)
library(ggplot2)
library(ggpubr)
library(paletteer)
library(rstatix)
library(tidyverse)
library(ggthemes)
library(ggprism)
library(ggsci)

#给热图添加样本注释信息：
group <- ifelse(str_sub(colnames(resm),14,15) < 10,"tumor","normal")

group <- as.data.frame(group_list)
colnames(group) <- "group" 
annotation <- data.frame(group)
row.names(annotation) <- colnames(resm)
table(annotation)
head(annotation)

pheatmap(resm,
         show_colnames = F,
         annotation_col = annotation,
         fontsize = 10)
#宽数据转换为长数据：
df <- resm %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(df)

#将value根据样本转换为百分比形式(新增一列)：
dff <- df %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dff)

#指定绘图顺序：
dff$cell_type <- factor(dff$cell_type,levels = unique(rownames(resm)))

#自定义主题：
mytheme <- theme(axis.title = element_text(size = 18),
                 axis.text.x = element_text(size = 15),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 18,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.position = "bottom")

#配色挑选
d_palettes <- palettes_d_names #查看离散型色板(paletteer包)
c_palettes <- palettes_c_names #连续型

p <- ggplot(dff,aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_paletteer_d("khroma::discreterainbow", direction = 1) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p
ggsave("cibersort免疫细胞差异.png", p, width = 15, height = 8, dpi = 600)

#指定箱线图排序（按相对丰度中位数从大到小）：
dff_arrange <- dff %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dff$cell_type <- factor(dff$cell_type,levels = unique(dff_arrange))

p1 <- ggplot(dff,aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_paletteer_d("khroma::discreterainbow", direction = 1) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p1

#新增分组列：
dff$group <- ifelse(str_sub(dff$sample,14,15)<10,"tumor","normal")

#新增分组
group <- cbind(colnames(resm),group)
colnames(group)[1] <- "sample"
dff <- merge(dff,group,by="sample",all = TRUE)

p2 <- ggplot(dff, aes(x = cell_type,y = proportion,fill = group)) +
  geom_boxplot(color = "black",alpha = 0.9,outlier.shape = 21,outlier.size = 1.2) +
  #scale_fill_manual(values = c("#4979b6","#d9352a")) +
  scale_fill_aaas()+
  labs(x = "cell type", y = "proportion") +
  theme_bw() + 
  mytheme + 
  theme_prism(border = T, base_rect_size = 1)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 8))+
  theme(legend.position = "top")
p2

#使用t test或wilcox test进行两两比较(T检验为例)：
t <- t_test(group_by(dff, cell_type), proportion ~ group)
t <- wilcox_test(group_by(dff, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr') #p值矫正；
tj

#根据p.adj添加显著性标记符号；
tj <- add_significance(tj, 'p.adj')
tj

#在图表中添加 p 值或者显著性标记；
lab <- add_xy_position(tj, x = 'cell_type')

p3 <- ggboxplot(dff, x = "cell_type", y = "proportion",
                fill = "group", color = "black",
                width=0.7, alpha = 0.9,
                outlier.shape = 21, outlier.size = 1.2) +
  #scale_fill_manual(values = c("#4979b6","#d9352a")) +
  labs(x = "cell type", y = "proportion") +
  scale_fill_aaas()+
  theme_bw() + 
  theme_prism(border = T, base_rect_size = 1)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 10)) +
  theme(legend.position = "top")+
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=5, bracket.size=1, tip.length = 0.01,
                     y.position = 0.75)
p3
ggsave("cibersort免疫细胞分组差异.png", p3, width = 12, height = 8, dpi = 600)

####小提琴图
#定义细胞顺序
dff$cell_type <- factor(dff$cell_type,levels = unique(rownames(resm)))
#取每种细胞类型的最大占比
dff_line <- dff %>%
  group_by(group,cell_type) %>%
  summarise("Mx" = max(proportion)) %>%
  group_by(cell_type) %>%
  summarise(lable = max(Mx))%>%
  arrange(cell_type)
#排序
dff_line <- dff_line %>% arrange(cell_type)
#定义小短线位置
x=c()
xend=c()
for(i in 1:length(unique(dff_line$cell_type))){
  x[i]=i-0.2
  xend[i]=i+0.2
}
y=dff_line$lable+0.01
yend=dff_line$lable+0.01
#整理成数据框
line <- data.frame(x,xend,y,yend)
#ggpubr算不了就手动计算P值
p_values <- dff %>%
  group_by(cell_type) %>%
  summarize(p_value = wilcox.test(proportion ~ group)$p.value) %>%
  mutate(p_label = ifelse(p_value < 0.001, sprintf("p < 0.001"), sprintf("p = %.3f", p_value)))

#画图
p_violin <- ggplot()+
  geom_violin(data=dff,aes(x=cell_type,y=proportion,fill=group),
              position = position_dodge(1),scale = "width")+
  geom_boxplot(data=dff,aes(x=cell_type,y=proportion,fill=group),
               position = position_dodge(1),width=0.4, outlier.shape = NA)+
  stat_summary(data=dff,aes(x=cell_type,y=proportion,fill=group),
               fun.y = "median",geom = "point",shape = 16,size = 2,color = "white",position = position_dodge(1))+
  scale_fill_aaas()+
  #映射出现报错，应该是是ggpubr跟ggplot有冲突了，看后面更新会不会修复
  #stat_compare_means(data=dff,aes(x=cell_type,y=proportion,group=group,label=ifelse(p<1.e-2,sprintf("p = %2.1e",as.numeric(..p.format..)),sprintf("p = %5.4f",as.numeric(..p.format..)))),method="wilcox.test",label.y=dff_line$lable+0.03,size=5)+
  geom_text(data=p_values,aes(x=cell_type,label = p_label), y = dff_line$lable + 0.03, 
            size = 2.5,fontface = "bold") +
  geom_segment(data=line,aes(x=x,y=y,xend=xend,yend=yend),linewidth=0.7)+
  labs(x = "", fill = "")+
  #theme_few(base_size = 12)+
  theme_prism(border = T, base_rect_size = 1,base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,colour = "black"),
        legend.position = c(0,1), legend.justification = c(0,1),legend.text = element_text(face = "bold"),
        legend.background = element_rect(fill = 'white', colour = 'black',linewidth=0.6))
p_violin
ggsave("cibersort-violin分组差异.pdf",p_violin,width = 25,height = 18,units = "cm",dpi = 600)

library(corrplot)
#计算相关性系数：
resmcor <- cor(t(resm))
#添加显著性标记:
resmorp <- cor.mtest(resmcor, conf.level = .95) #使用cor.mtest做显著性检验;
#提取p值矩阵；
p.mat <- resmorp$p

png(filename = "cibersort-corrplot.png",width = 2600,height = 2600,res=300)
corrplot(resmcor, method = "color",order = "hclust",
         type="upper", tl.pos = "lt",tl.col="black",tl.cex = 0.8,
         p.mat = resmorp$p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.8, pch.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(50))
corrplot(resmcor,method = "color",order = "hclust",
         type="lower",add = T,
         tl.pos="n",cl.pos="n",diag = F,
         col = colorRampPalette(c("blue", "white", "red"))(50),
         number.cex = 0.7,addCoef.col = "black")
dev.off()


#选取目标基因
load("C:/Users/11146/Desktop/R包/TCGA/int_gene.Rdata")
load("KICH_exp_tpm.Rdata")
target_genes <- c("MX1","IFIT3","GBP1","SERPING1","IFIT1","IFI44","ISG15","IFI27")
target_genes <- "IFI27"
#对TPM取log：
exp_log2 <- log2(KICH_exp_tpm + 1)
genes <- int_p1
exp_genes <- exp_log2[genes,]

load("exp.Rdata")
#提取目标基因矩阵
exp_genes <- exp[target_genes,]
#去除缺失值的行
exp_genes <- na.omit(exp_genes)
#将免疫细胞丰度矩阵和目的基因表达量矩阵合并：
rb <- rbind(resm,exp_genes)
#重新计算相关性系数：
rbcor <- cor(t(rb))
#计算显著性差异：
rbcorp <- cor.mtest(rbcor, conf.level = .95)
#提取p值矩阵：
p.mat2 <- rbcorp$p

#切割相关性矩阵：
split <- rbcor[1:nrow(resm), #行取免疫细胞所在的行
               (nrow(rb)-length(target_genes)+1):nrow(rb)] #列取目的基因所在的行
View(split) #行名免疫细胞，列名目的基因；
write.csv(split,"DEGs与免疫细胞相关性.csv", row.names = T, quote = F)

#切割p值矩阵:
splitp <- p.mat2[1:nrow(resm), #行取免疫细胞所在的行
                 (nrow(rb)-length(target_genes)+1):nrow(rb)] #列取目的基因所在的行
View(splitp)

mark <- matrix(ifelse(splitp < 0.001, "***", 
                      ifelse(splitp < 0.01, "**", 
                             ifelse(splitp < 0.05, "*", ""))),
               nrow = nrow(splitp))

save(split,splitp,mark,file = c("cibersort-split-splitp.Rdata"))

load("cibersort-split-splitp.Rdata")

col2 <- colorRampPalette(c("blue","white", "red"))(95)
pheatmap(split,
         display_numbers = mark,
         number_color = "black",
         fontsize_number = 8,
         cellwidth = 20,
         cellheight = 15,
         color = col2,
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_col = 8,
         fontsize_row = 8,
         angle_col = 45,
         filename = "cibersort免疫细胞-基因相关性热图.pdf",
         width = 20,height = 15)


library(dplyr)
library(tidyr)
library(tidyverse)
# 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Sample”。
#并赋值给cibersort_raw。
cibersort_raw <- read.table("CIBERSORT-Results.txt", header = T, sep = '\t') %>% 
  rename("Sample" = "Mixture") %>%
  select(-c("P.value","Correlation","RMSE"))
# 将cibersort_raw第一列变为列名后赋值给cibersort_tidy。
cibersort_tidy <- cibersort_raw %>%
  remove_rownames() %>%
  column_to_rownames("Sample")
# 筛选出0值太多的一些细胞。
flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) < 
                dim(cibersort_tidy)[1]/2)
# 留下在大部分样本中有所表达的细胞。
cibersort_tidy <- cibersort_raw[,which(flag)] %>%
  as.matrix() %>%
  t()
# breaks用来定义数值和颜色的对应关系。
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))

# 将CIBERSORT_Result进行可视化
#使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值
df <- cibersort_raw %>%
  gather(key = Cell_type,value = Proportion,2:23)

df_tumor <- df[substr(df$Sample, 14, 15) < 10, ]
df_normal <- df[substr(df$Sample, 14, 15) > 10, ]

library(RColorBrewer)
library(ggplot2)
mypalette <- colorRampPalette(brewer.pal(5,"Set1"))
#绘制柱状图
p1 <- ggplot(df_tumor,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell type",x = "tumor",y = "Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values=mypalette(23))+
  theme(legend.key.size = unit(0.3,'cm'))
p1
ggsave("cibersort-tumor柱状图.png",p1,width = 20,height = 14,units = "cm",dpi = 600)

p2 <- ggplot(df_normal,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell type",x = "normal",y = "Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values=mypalette(23))+
  theme(legend.key.size = unit(0.3,'cm'))
p2
ggsave("cibersort-normal柱状图.png",p2,width = 15,height = 14,units = "cm",dpi = 600)


#整合在一起
df$group <- ifelse(str_sub(df$Sample,14,15) < 10,"Tumor","Normal")

group <- as.data.frame(group_list)
colnames(group) <- "group" 
group <- cbind(row.names(res),group)
colnames(group)[1] <- "Sample"
df <- merge(df,group,by="Sample",all = TRUE)

#按照分组进行排列
df <- df %>% mutate(Sample = factor(Sample, levels = unique(df$Sample[order(df$group)])))

p <- ggplot(df,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell type",x = "Sample",y = "Proportion") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values=mypalette(23))+
  theme(legend.key.size = unit(0.3,'cm'))
p
ggsave("cibersort-sample柱状图.png",p,width = 20,height = 14,units = "cm",dpi = 600)

#定义分面颜色
library(ggh4x)
ridiculous_strips <- strip_themed(background_x = elem_list_rect(fill = c("#3b4992","#ee0000")))
#或者用分面来做
p <- ggplot(df, aes(x = Sample, y = Proportion, fill = Cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.7) +
  labs(fill = "Cell type", x = "", y = "Proportion") +
  theme_prism(base_size = 10) +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = mypalette(23)) +
  theme(legend.key.size = unit(0.3, 'cm')) +
  #facet_grid(. ~ group, scales = "free_x", space = "free_x", switch = "x")+
  facet_grid2(.~group, scales = "free_x", space = "free_x", switch = "x",strip = ridiculous_strips)+
  theme(strip.text = element_text(face = "bold",colour = "white"))+
  theme(panel.spacing.x = unit(0, "cm"))+#调整分面间距
  theme(legend.text = element_text(face = "bold"))+
  guides(fill = guide_legend(ncol = 1))
p
ggsave("cibersort-sample柱状图分面.pdf",p,width = 25,height = 15,units = "cm",dpi = 600)


#PCA分析
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
load("exp.Rdata")
cibersort_raw <- read.table("CIBERSORT-Results.txt", row.names = 1, header = T, sep = '\t') %>%
  select(-c("P.value","Correlation","RMSE"))
cibersort_raw <- t(cibersort_raw)

#group_list <- ifelse(as.numeric(str_sub(colnames(cibersort_raw),14,15)) < 10,'tumor','normal')
#group_list <- factor(group_list,levels = c("normal","tumor"))

#if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)
draw_pca(cibersort_raw,group_list)

#美化
cibersort_raw <- read.table("CIBERSORT-Results.txt", row.names = 1, header = T, sep = '\t') %>%
  select(-c("P.value","Correlation","RMSE"))
#group <- ifelse(as.numeric(str_sub(colnames(cibersort_raw),14,15)) < 10,'tumor','normal')
group <- as.data.frame(group_list)
colnames(group) <- "group" 
group <- cbind(row.names(cibersort_raw),group)
colnames(group)[1] <- "Sample"
group$group <- factor(group$group,levels = c("control","HIV"))

# 主成分分析
pca <- prcomp(cibersort_raw,center = T,scale. = F)
#center 一个逻辑值，控制变量是否应该移位到零中心
#scale. 一个逻辑值，控制是否对数据进行标准化
#计算pc1、pc2
a <- summary(pca)
tmp <- a$importance
pc1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
pc2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
# 获取距离矩阵
pca_mat <- data.frame(pc1=pca$x[,1],pc2=pca$x[,2],sample=group$Sample,group=group$group)


library(ggplot2)
p <- ggplot(pca_mat,aes(x=pc1,y=pc2,colour=group)) +
  geom_point(size=2) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
  labs(x=paste("PC1(",pc1,"%)",sep=""),y=paste("PC2(",pc2,"%)",sep=""),
       title = "") +
  scale_fill_manual(values=c("#3C5488B2", "#DC0000B2"))+
  scale_color_manual(values=c("#3C5488B2", "#DC0000B2"))+
  geom_hline(yintercept=0,linetype='dotdash',linewidth=0.8,color='grey') +
  geom_vline(xintercept=0,linetype='dotdash',linewidth=0.8,color='grey') +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5))
p
ggsave("cibersort-PCA-ggplot2.png",p,width = 15,height = 12,units = "cm",dpi = 600)

