#GSEA和GSVA都是基于对基因的某一个值的排序来进行富集分析
#GSEA主要是用case和control之间的差异倍数或信噪比来进行排序（一次处理两个样本）
#GSVA则不需要做对比，而是对每个样本或单个细胞按基因的表达量进行单独排序
#然后将富集分数的值做个标准化，再在不同的样本间对比
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
rm(list = ls())
options(stringsAsFactors = F)#否则所有的字符在数据框中都会被R默认设置为factor

setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE33463")
load('exp.Rdata') #包含 tpm counts group_list gl
#创建新文件夹
dir.create("./GSVA")
setwd("./GSVA")

d <- "C:/Users/maihuanzhuo/Desktop/R包/MSigDB"
## 对 MigDB( Molecular Signatures Database)中的基因集做GSVA分析
gmtfs <- list.files(d,pattern = 'symbols.gmt')  # 路径下所有结尾为symbols.gmt文件
gmtfs
hallmark_list <- getGmt(file.path(d,gmtfs[31]))#"GSEABase"
kegg_list <- getGmt(file.path(d,gmtfs[5])) 
go_list <- getGmt(file.path(d,gmtfs[24])) 

## msigdbr包提取下载 一般尝试KEGG和GO做GSVA分析
##KEGG
KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
                        category = "C2",
                        subcategory = "CP:KEGG") 
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组

##GO
GO_df_all <- msigdbr(species = "Homo sapiens",
                     category = "C5")  
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##按照gs_name给gene_symbol分组

##hallmark
hallmark_df_all <- msigdbr(species = "Homo sapiens",
                           category = "H")
hallmark_df <- dplyr::select(hallmark_df_all,gs_name,gs_exact_source,gene_symbol)
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name) ##按照gs_name给gene_symbol分组

#GSVA
#GSVA算法需要处理logCPM, logRPKM,logTPM数据或counts数据的矩阵####
#GSVA的矩阵是经过log2标准化处理的芯片数据或者是rna-seq的count data
#dat <- as.matrix(counts)
#dat <- as.matrix(log2(edgeR::cpm(counts))+1)
#dat <- as.matrix(log2(tpm+1))
str(exp)
class(exp)

gsva_mat <- gsva(expr = exp, 
                 gset.idx.list = hallmark_list, #kegg_list
                 method = "gsva",
                 kcdf = "Gaussian" ,#Gaussian (针对芯片表达矩阵，log CPMs, log RPKMs, log-TPM标准化后RNA-seq的表达矩阵)，
                 #可选Poisson (针对未标准化的RNA-seq的count表达矩阵) ;
                 verbose = T, 
                 #min.sz = 5, max.sz = 500,
                 parallel.sz = parallel::detectCores())#调用所有核
#每个样品都会有对应通路的GSVA评分
write.csv(gsva_mat,"gsva_hallmark_matrix.csv")

#目标基因与
target_genes <- c("MX1","IFIT3","GBP1","SERPING1","IFIT1","IFI44","ISG15","IFI27")
# 从基因表达矩阵中提取目标基因的表达矩阵
target_gene_matrix <- exp[target_genes, ]

#limma差异分析
#得到GSVA评分的矩阵后，我们需要利用limma包进行pathway通路的差异分析
#但不需要进行 limma-trend 或 voom的步骤
#### 进行limma差异处理 ####
##设定 实验组exp / 对照组ctr
group_list
experiment <- "PAH"
ctr <- "control"

design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(gsva_mat)
contrast.matrix <- makeContrasts(contrasts = paste0(experiment,'-',ctr),  #"exp/ctrl"
                                 levels = design)

fit1 <- lmFit(gsva_mat,design)#拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)#贝叶斯修正

summary(decideTests(efit, p.value = 0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef = paste0(experiment,'-',ctr), n = Inf)
degs <- na.omit(tempOutput) 
write.csv(degs,"gsva_hallmark_degs.results.csv")


#可视化
##设置筛选阈值
padj_cutoff = 0.05
log2FC_cutoff = log2(1)#1

keep <- rownames(degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC) > log2FC_cutoff, ])
length(keep)
dat <- gsva_mat[keep[1:28],] #选取前50进行展示
#热图
pheatmap::pheatmap(dat, 
                   fontsize_row = 8,
                   height = 10,
                   width = 18,
                   #annotation_col = group_list,
                   show_colnames = F,
                   show_rownames = T,
                   filename = paste0('GSVA_hallmark_heatmap.pdf'))

#火山图
degs$significance  <- as.factor(ifelse(degs$adj.P.Val < padj_cutoff & abs(degs$logFC) > log2FC_cutoff,
                                       ifelse(degs$logFC > log2FC_cutoff ,'UP','DOWN'),'NOT'))

this_title <- paste0(' Up :  ',nrow(degs[degs$significance =='UP',]) ,
                     '\n Down : ',nrow(degs[degs$significance =='DOWN',]),
                     '\n adj.P.Val <= ',padj_cutoff,
                     '\n FoldChange >= ',round(2^log2FC_cutoff,3))

g <- ggplot(data=degs, 
            aes(x=logFC, y=-log10(adj.P.Val),
                color=significance)) +
  #点和背景
  geom_point(alpha=0.4, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( adj.P.Val )") +
  #标题文本
  ggtitle( this_title ) +
  #分区颜色                   
  scale_colour_manual(values = c('blue','grey','red'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(),
        legend.position="right")
g
ggsave(g,filename = 'GSVA_go_volcano_padj.pdf',width = 8,height = 7.5, dpi = 600)


#### 发散条形图绘制 ####
library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(ggthemes)
library(ggprism)
p_cutoff=0.05
degs <- gsva_kegg_degs  #载入gsva的差异分析结果
Diff <- rbind(subset(degs,logFC>0)[1:20,], subset(degs,logFC<0)[1:8,]) #选择上下调前20通路     
dat_plot <- data.frame(id = row.names(Diff),
                       p = Diff$P.Value,
                       lgfc = Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$id <- str_replace(dat_plot$id, "KEGG_","");dat_plot$id[1:10]
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels = c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## 设置不同标签数量
low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
high1 <- nrow(dat_plot)

p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, fill = threshold)) +
  geom_col()+
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#36638a','Not'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('-log10(P.Value) of GSVA score') + 
  guides(fill="none")+
  theme_prism(border = T) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black') + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 黑色标签
p
ggsave("GSVA_barplot_pvalue.pdf",p,width = 15,height  = 15,dpi = 600)
