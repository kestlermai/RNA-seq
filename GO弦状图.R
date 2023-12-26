#install.packages("GOplot")
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOplot)
library(RColorBrewer)
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE77939")
load("nrDEG_limma_voom.Rdata")
#logFC=0.585 or 1 or 1.5 or 2    adj.P还是P.Value
nrDEG_limma_voom_signif <- nrDEG_limma_voom %>% filter(abs(logFC) > 0.585) %>% filter(P.Value < 0.05)
write.csv(nrDEG_limma_voom_signif,"limma-signif.csv")
nrDEG_limma_voom_signif <- cbind(gene = rownames(nrDEG_limma_voom_signif), nrDEG_limma_voom_signif)
write.table(nrDEG_limma_voom_signif[,c(1,2,5)],"GSE77939_limma_signif.txt",sep="\t",quote=F,row.names = F)
#制作成对应的格式
genelist <- nrDEG_limma_voom_signif
row.names(genelist) <- 1:nrow(genelist)
colnames(genelist)[1] <- "ID"
#转换ENTREZID
gene_id <- bitr(genelist$gene, 
                fromType = "SYMBOL", #需要转换的类型
                toType = c("ENTREZID"), #需要转换为的类型
                OrgDb = org.Hs.eg.db) #注释包

GO_all <- enrichGO(gene = gene_id$ENTREZID, 
                        OrgDb = org.Hs.eg.db, 
                        ont = "all", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = T) 
GO_all_result <- GO_all@result
#制作成对应的格式
row.names(GO_all_result) <- 1:nrow(GO_all_result)
colnames(GO_all_result)[1] <- "Category"
colnames(GO_all_result)[3] <- "Term"
colnames(GO_all_result)[7] <- "adj_pval"
colnames(GO_all_result)[9] <- "Genes"
GO_all_result$Genes <- gsub("/", ",", GO_all_result$Genes)
GO_all_result <- GO_all_result[,c(1,2,3,9,7)]

KEGG_diff <- enrichKEGG(gene = gene_id$ENTREZID,
                        organism = "hsa",#物种，Homo sapiens (human)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
KEGG_result <- KEGG_diff@result

#将两个表关联起来
circ <- circle_dat(GO_all_result, genelist)
#circ可以直接画circle圆圈图
GOCircle(circ,rad1=2, #内环半径
         rad2=3, #外环半径
         label.size= 5, #标签大小
         label.fontface = 'bold', #标签字体
         nsub=10, #显示的terms数，前10个。（画图前需要先把数据整理好，想要哪些term）
         zsc.col = c('red', 'white', "green"), # z-score颜色
         lfc.col = c('red', 'green')) # 基因up或down的颜色
##Chord弦图
#先创建二进制矩阵，有点类似热图数据
chord <- chord_dat(circ, 
                   genelist[1:100,], #选择展示的基因
                   GO_all_result$Term[1:10])#选择展示的term
head(chord)
#绘图
pdf("chord.pdf",height = 12,width = 12)
GOChord(chord, space = 0.02, #弦的间隔
        gene.order = 'logFC', #可以是-log（Pvalue）值
        gene.space = 0.25, gene.size = 5,
        border.size = NULL, #彩虹边框大小
        lfc.col = c('red','white','blue'), #自定义logFC 颜色
        ribbon.col = brewer.pal(10,"Paired"))#更改配色
dev.off()

#去除adj_P大于0.05
genelist <- genelist_test %>%tibble::rownames_to_column(var="ID") 
df <- as.data.frame(genelist_test)
df$Regulate =ifelse(df$logFC <=-1& df$padj_DESeq2<=0.05,'Down',
                    ifelse(df$logFC >=1 & df$padj_DESeq2<=0.05,'Up','Stream'))
new_df <- df[which(df$Regulate=='Up'|df$Regulate=='Down'),]
Genelist <- genelist[genelist$ID %in% new_df$name,]

#
library(org.Hs.eg.db)
library(clusterProfiler)
gene_symbol <- c("IFI27","IFI6","IFI27L1","IFI27L2","NR4A1","IRF9","STAT2","STAT1","ISG15","IFITM1","MX1",
                 "OAS1","IFI44L","IFIT3","IFIT1","IFI35","BST2","RSAD2","OAS2","XAF1","MX2")
gene_symbol <- data.frame(gene_symbol)
gene_id <- bitr(gene_symbol$gene_symbol, 
                fromType = "SYMBOL", #需要转换的类型
                toType = c("ENTREZID"), #需要转换为的类型
                OrgDb = org.Hs.eg.db) #注释包
head(gene_id)

GO_all_diff <- enrichGO(gene = gene_id$ENTREZID, 
                        OrgDb = org.Hs.eg.db, 
                        ont = "all", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = T) 
GO_all_result <- GO_all_diff@result


