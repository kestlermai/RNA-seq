#DEGs差异分析
#先画PCA看看能不能做差异分析
setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE77939")

library(limma)
library(tinyarray)
library(dplyr)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(ggplot2)

load("exp.Rdata")
dim(exp)
dat <- as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = F)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, 
                         palette = c("#0072B5FF","#BC3C29FF"),
                         addEllipses = T, # Concentration ellipses
                         legend.title = "Group")
pca_plot
ggsave("pca_plot.pdf",pca_plot,width = 15,height = 12,units = "cm",dpi = 600)

#PCA分析
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE77939")
load("exp.Rdata")
group <- as.data.frame(group_list)
colnames(group) <- "group" 
exp <- t(exp)
group <- cbind(row.names(exp),group)
colnames(group)[1] <- "Sample"
group$group <- factor(group$group,levels = c("control","HIV"))

# 主成分分析
pca <- prcomp(exp,center = T,scale. = T)
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
ggsave("Sample-PCA-ggplot2.tiff",p,width = 15,height = 12,units = "cm",dpi = 600)

#绘制热图
cg <- names(tail(sort(apply(exp,1,sd)),1000))
n <- exp[cg,]
annotation_col <- data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 

t <- t(scale(t(n)))
table(abs(t)>2)
t[t>=2]=2
t[t<=-2]=-2

pheatmap(n,
         border_color = 'black',
         show_colnames = F,
         show_rownames = F,
         cluster_col= T,
         cluster_rows= T,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100))


setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE703")
load("exp.Rdata")
#limma差异分析比较
library(limma)
library(dplyr)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp)
#Cancer写左，Normal写右
contrast.matrix<-makeContrasts(PAH-control,levels = design)
contrast.matrix
#非线性最小二乘法
fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
#用经验贝叶斯调整t-test中方差的部分
fit2 <- eBayes(fit2)
#voom
nrDEG_limma_voom <- topTable(fit2, coef=1, n=Inf)
nrDEG_limma_voom <- na.omit(nrDEG_limma_voom)
write.csv(nrDEG_limma_voom,"limma.csv")
save(nrDEG_limma_voom,file = 'nrDEG_limma_voom.Rdata')
#筛选显著性差异的基因
setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE77939")
load("nrDEG_limma_voom.Rdata")
#logFC=0.585 or 1 or 1.5 or 2    adj.P还是P.Value
nrDEG_limma_voom_signif <- nrDEG_limma_voom %>% filter(abs(logFC) > 0.585) %>% filter(P.Value < 0.05)
write.csv(nrDEG_limma_voom_signif,"limma-signif_0.585.csv")
nrDEG_limma_voom_signif <- cbind(ID = rownames(nrDEG_limma_voom_signif), nrDEG_limma_voom_signif)
write.table(nrDEG_limma_voom_signif[,c(1,2,5)],"GSE77939_limma_signif.txt",sep="\t",quote=F,row.names = F)
#进行GO/KEGG富集分析
genelist <- nrDEG_limma_voom_signif
row.names(genelist) <- 1:nrow(genelist)

## 火山图绘制R包
#devtools::install_github("BioSenior/ggVolcano")
nrDEG_limma_voom$gene <- row.names(nrDEG_limma_voom)
library(ggVolcano)
data <- add_regulate(nrDEG_limma_voom, 
                     log2FC_name = "logFC",
                     fdr_name = "P.Value",
                     log2FC = 0.585, 
                     fdr = 0.05)
ggvolcano(data, 
          x = "log2FoldChange",
          y = "padj",
          label = "gene", 
          fills = c("#e94234","grey","#269846"),
          colors = c("#e94234","grey","#269846"),
          label_number = 0, # 调整显示的标签数量
          log2FC_cut = 0.585,
          FDR_cut = 0.05,
          output = FALSE,
          legend_position="DR")
ggsave(file = "volcano.pdf", width = 10, height = 8)

#DEGs热图
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE703")
load("nrDEG_limma_voom.Rdata")
load("exp.Rdata")
nrDEG_limma_voom_signif <- nrDEG_limma_voom %>% filter(abs(logFC) > 0.585) %>% filter(P.Value < 0.05)
DEG_gene_expr <- exp[rownames(nrDEG_limma_voom_signif),]

annotation_col <- data.frame(group=group_list)
rownames(annotation_col)=colnames(DEG_gene_expr) 

pheatmap(DEG_gene_expr,
         cluster_row= T,
         cluster_col= T,
         #annotation_col=annotation_col,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         breaks = seq(-3.5,3.5,length.out = 100),
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         fontsize = 10, #文字大小
         show_rownames = F,
         filename = "diff_heatmap.pdf",
         width = 12,height = 8)




#RRA算法对基因进行排序
#BiocManager::install("RobustRankAggreg")
library(RobustRankAggreg)
setwd("C:/Users/11146/Desktop/双疾病生信数据/RRA-PAH")
files <- c("GSE703_limma_signif.txt","GSE19617_limma_signif.txt","GSE33877_limma_signif.txt",
           "GSE77939_limma_signif.txt","GSE131793_limma_signif.txt")

files <- c("GSE703_limma_signif.txt","GSE19617_limma_signif.txt")

library(readr)
library(dplyr)
upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile,header=T,sep = '\t',quote = '') # 注意文件读取
  header=unlist(strsplit(inputFile,"_"))
  downList[[header[1]]]=as.vector(rt[,1])
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]#只留下了rt的前两列，保留下的是基因名和LogFC
  colnames(fcCol)=c("gene",header[[1]])#第一列为行名
  allFCList[[header[1]]]=fcCol
}

padj=0.05
logFC=1

#合并所有差异基因logFC
mergeLe=function(x,y){
  merge(x,y,by="gene",all=T)}
newTab=Reduce(mergeLe,allFCList)#两个表格合并
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0

#筛选共同up基因
upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="bonferroni")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
write.table(upXls,file="up.xls",sep="\t",quote=F,row.names=F)
upSig=upXls[(upXls$Pvalue <padj & upXls$logFC>logFC),]#在这一步把adjPvalue改成Pvalue
write.table(upSig,file="upSig.xls",sep="\t",quote=F,row.names=F)

#筛选共同下调基
downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))
write.table(downXls,file="down.xls",sep="\t",quote=F,row.names=F)
downSig=downXls[(downXls$Pvalue<padj & downXls$logFC< -logFC),]
write.table(downSig,file="downSig.xls",sep="\t",quote=F,row.names=F)

#合并共同上下调基
allSig = rbind(upSig,downSig)
colnames(allSig)
allSig = allSig[,c("Name","logFC")]
write.table(allSig,file = 'allSign.xls',sep = '\t',quote = F)
#logFC.tiff
hminput=newTab[c(as.vector(upSig[1:10,1]),as.vector(downSig[1:10,1])),]
hminput=newTab[c(as.vector(upXls[1:10,1]),as.vector(downXls[1:10,1])),]
colnames(hminput) <- c("GSE77939", "GSE33877","GSE703","GSE19617","GSE131793")
colnames(hminput) <- c("GSE703","GSE19617")
library(pheatmap)
tiff(file="logFC.tiff",width = 15,height = 20,units ="cm",compression="lzw",bg="white",res=400)
pheatmap(hminput,display_numbers = TRUE,
         fontsize_row=10,
         fontsize_col=12,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = FALSE,cluster_rows = FALSE, )
dev.off()


