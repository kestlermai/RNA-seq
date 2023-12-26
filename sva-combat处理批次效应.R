#BiocManager::install("sva")
#SVA包处理批次效应前并不需要limma包校正
library(sva)
#PAH
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
load("exp.Rdata")
exp1 <- exp
group1 <- group_list
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE703")
load("exp.Rdata")
exp2 <- exp
exp2 <- as.matrix(exp2)
group2 <- group_list
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE131793")
load("exp.Rdata")
exp3 <- exp
group3 <- group_list

sameSample <- intersect(rownames(exp1), rownames(exp2))
sameSample <- as.data.frame(sameSample)
sameSample <- intersect(rownames(exp3), sameSample$sameSample)
gene_exp1 <- exp1[sameSample,,drop=F]
gene_exp2 <- exp2[sameSample,,drop=F]
gene_exp3 <- exp3[sameSample,,drop=F]

bindgeo <- cbind(gene_exp1,gene_exp2,gene_exp3)


group1 <- as.matrix(group1)
colnames(group1) <- "group"
group2 <- as.matrix(group2)
colnames(group2) <- "group"
group3 <- as.matrix(group3)
colnames(group3) <- "group"

talgroup <- as.data.frame(rbind(group1,group2,group3))
talgroup_list <- factor(talgroup$group,levels = c("control","PAH"))
setwd("C:/Users/11146/Desktop/双疾病生信数据/cbind-PAH")
save(bindgeo,talgroup,file = "exp.Rdata")

#箱式图
png("boxplot.png", width=1200, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
boxplot(bindgeo, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()
#绘制去批次前PCA图
library(ggplot2)
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(as.data.frame(t(bindgeo)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = talgroup_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="prenormal_PCA.png")

batchType <- c(rep(1,113),rep(2,20),rep(3,20))
modType <- c(rep("control",41),rep("PAH",72),rep("PAH",14),rep("contorl",6),
          rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),
          rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),
          rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),rep("PAH",1),rep("control",1),
          rep("PAH",1),rep("control",1))
mod <- model.matrix(~modType)
library(sva)
bindgeo_ComBat <- ComBat(dat=bindgeo, batch=batchType, #使用ComBat法去批次
                      mod=mod, par.prior=T)
#box
png("boxplot_combat.png", width=1200, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
boxplot(bindgeo_ComBat, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()
#PCA
dat.pca2 <- PCA(as.data.frame(t(bindgeo_ComBat)), graph = FALSE)
pca_plot2 <- fviz_pca_ind(dat.pca2,
                          geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                          col.ind = talgroup_list,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot2
ggsave(plot = pca_plot2,filename ="afternormal_PCA.pdf")

#热图
cg <- names(tail(sort(apply(bindgeo_ComBat,1,sd)),250))
n <- bindgeo_ComBat[cg,]
annotation_col <- data.frame(group=talgroup_list)
rownames(annotation_col)=colnames(n) 
p <- pheatmap(n,
             border_color = 'black',
             show_colnames = F,
             show_rownames = F,
             cluster_col= T,
             cluster_rows= T,
             annotation_col=annotation_col,
             scale = "row",
             breaks = seq(-3,3,length.out = 100))
ggsave(plot = p,filename ="afternormal_heatmap.pdf")

save(bindgeo_ComBat,talgroup,talgroup_list,file = "exp-combat.Rdata")

#HIV
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE140713")
load("exp.Rdata")
exp1 <- exp
group1 <- group_list
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33877")
load("exp.Rdata")
exp2 <- exp
group2 <- group_list
sameSample <- intersect(rownames(exp1), rownames(exp2))
gene_exp1 <- exp1[sameSample,,drop=F]
gene_exp2 <- exp2[sameSample,,drop=F]
bindgeo <- cbind(gene_exp1,gene_exp2)
group1 <- as.matrix(group1)
colnames(group1) <- "group"
group2 <- as.matrix(group2)
colnames(group2) <- "group"
talgroup <- as.data.frame(rbind(group1,group2))
talgroup_list <- factor(talgroup$group,levels = c("control","HIV"))
setwd("C:/Users/11146/Desktop/双疾病生信数据/cbind-HIV")
save(bindgeo,talgroup,file = "exp.Rdata")

png("boxplot.png", width=1200, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
boxplot(bindgeo, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()
#绘制去批次前PCA图
library(ggplot2)
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(as.data.frame(t(bindgeo)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = talgroup_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="prenormal_PCA.pdf")

batchType <- c(rep(1,57),rep(2,12))
modType <- c(rep("HIV",50),rep("control",7),
             rep("contorl",1),rep("HIV",1),
             rep("contorl",1),rep("HIV",1),
             rep("contorl",1),rep("HIV",1),
             rep("contorl",1),rep("HIV",1),
             rep("contorl",1),rep("HIV",1),
             rep("contorl",1),rep("HIV",1)
             )
mod <- model.matrix(~modType)
library(sva)
bindgeo_ComBat <- ComBat(dat=bindgeo, batch=batchType, #使用ComBat法去批次
                         mod=mod, par.prior=T)
#box
png("boxplot_combat.png", width=1200, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
boxplot(bindgeo_ComBat, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()

#PCA
dat.pca2 <- PCA(as.data.frame(t(bindgeo_ComBat)), graph = FALSE)
pca_plot2 <- fviz_pca_ind(dat.pca2,
                          geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                          col.ind = talgroup_list,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot2
ggsave(plot = pca_plot2,filename ="afternormal_PCA.pdf")
#热图
cg <- names(tail(sort(apply(bindgeo_ComBat,1,sd)),800))
n <- bindgeo_ComBat[cg,]
annotation_col <- data.frame(group=talgroup_list)
rownames(annotation_col)=colnames(n) 
p <- pheatmap(n,
         border_color = 'black',
         show_colnames = F,
         show_rownames = F,
         cluster_col= T,
         cluster_rows= T,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100))
ggsave(plot = p,filename ="afternormal_heatmap.pdf")

save(bindgeo_ComBat,talgroup,talgroup_list,file = "exp-combat.Rdata")

library(limma)
bindgeo_after <- bindgeo_ComBat
design=model.matrix(~talgroup_list)
fit=lmFit(bindgeo_after,design)#这里要注意，分组的样本与矩阵样本是否相符，不相符则去文件中调整然后读入
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
logFC=1
adj.P.Val = 0.05
k1 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="upanddown.csv")
