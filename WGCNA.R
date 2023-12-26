#WGCNA 
library(WGCNA)
library(doParallel)
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
load("exp.Rdata")
#（1）导入数据
#输入数据的格式是列为样本，行为基因的矩阵格式，所以不用置换矩阵
## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置。
#根据其方差选择了前25%的变异较大的基因，包括约5，000个基因进行WGCNA分析
exp_mt = t(exp[order(apply(exp,1,mad), decreasing = T)[1:5000],])

##(2)判断数据质量
gsg <- goodSamplesGenes(exp_mt,verbose = 3)
gsg$allOK
if (!gsg$allOK) {                    # 如果存在异常样本或基因
  if (sum(!gsg$goodGenes)>0)       # 异常的基因
    printFlush(paste("Removing genes:", 
                     paste(names(exp_mt)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)     # 异常的样本
    printFlush(paste("Removing samples:", 
                     paste(rownames(exp_mt)[!gsg$goodSamples], collapse = ",")));
  # 删除异常的样本和基因
  exp_mt = exp_mt[gsg$goodSamples, gsg$goodGenes]
}

##(3)判断数据质量--离群点样本
sampleTree = hclust(dist(exp_mt), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 70, col = "red") #根据实际情况而定
##如下图所示，存在一个显著离群点；剔除掉
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
table(clust)
# clust
# 0   1 
# 1 112
keepSamples = (clust==1)
exp_mt_f = exp_mt[keepSamples, ]
exp_mt_f <- exp_mt
row <- row.names(exp_mt_f)
row <- data.frame(row)
write.csv(row,"rowgroup.csv")
#分组数据
group <- read.table("group.txt",row.names = 1,header = T, sep = "\t")

sameSample=intersect(rownames(exp_mt_f), rownames(group))
exp_mt_f=exp_mt_f[sameSample,]
datTraits=group[sameSample,]
traitColors = numbers2colors(datTraits, signed = FALSE)

sampleTree2 = hclust(dist(exp_mt_f), method = "average")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree2, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

#选择合适的软阈值β
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(exp_mt_f, powerVector = powers, verbose = 5)
##结果可视化：
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.85#一般是0.9，可以降低为0.8，最低是0.8
### （1）是否符合幂律分布；
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
#横坐标为软阈值的梯度，第一幅图的纵坐标为无标度网络适应系数，越大越好；第二幅图的纵坐标为节点的平均连通度，越小越好。
#基因分布符合无尺度网络
#查看系统给我们推荐的软阈值
sft$powerEstimate
#4
#6

####----构建共表达网络----####
net <- blockwiseModules(
  exp_mt_f,
  power = 6,
  minModuleSize = 20,#模块的最少基因数
  corType = "pearson",#默认皮尔森
  networkType="unsigned",#计算邻接矩阵时，是否考虑正负相关性，默认unsigned
  TOMType = "unsigned",#计算TOM矩阵时，是否考虑正负相关性，默认unsigned
  mergeCutHeight = 0.25,#合并模块的阈值
  numericLabels = TRUE,#模块名是否为数字；若设置FALSE，表示映射为颜色名
  saveTOMs = TRUE,
  saveTOMFileBase = "PAH-TOM",
  verbose = 5
)

names(net)

#聚类树状图，color注释模块
#灰色部分太多则说明由于样本中基因共表达趋势不明显，可能需要调整基因过滤的方法
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# (1)计算模块特征值
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

nGenes = ncol(exp_mt_f)
nSamples = nrow(exp_mt_f)

# 用color labels重新计算MEs（Module Eigengenes:模块的第一主成分）
MEs0 = moduleEigengenes(exp_mt_f, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# (2)计算module与表型的相关性以及对应的P值
moduleTraitCor = cor(MEs, group, use = "p")
moduleTraitCor[1:2,1:2]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs))

#(3)可视化相关性与P值
sizeGrWindow(10,6)
# #设置热图上的文字（两行数字：第一行是模块与各种表型的相关系数；第二行是p值，signif 取有效数字
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 然后对moduleTraitCor画热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(group),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#相似模块进行网络构建
MEList = moduleEigengenes(exp_mt_f, colors = moduleColors)#【提取相似模块】
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#筛选hub基因
#计算基因与模块的相关性
nSamples <- dim(exp_mt_f)[1]
geneModuleMembership <- cor(exp_mt_f, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
#计算基因与表型的相关性
geneSignificanceCor <- cor(exp_mt_f, group, use = "p")
geneSignificanceP <- corPvalueStudent(geneSignificanceCor, nSamples)
#模块内基因与模块和表型之间的相关性
module <- "blue"
column <- paste0("ME", module)
moduleGenes <- names(net$colors)[which(moduleColors == module)]
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneSignificanceCor[moduleGenes, 1])
verboseScatterplot(
  MM, GS,
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for proliferating",
  main = paste("Module membership vs. gene significance\n"),
  abline = TRUE,
  pch = 21,
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = "black",
  bg = module
)

#推荐选取散点图中右上角部分的基因，即两个相关性均较大的基因。
#brown
moduleGenes[(GS > 0.4 & MM > 0.85)]
# "SCYB10";"BAL";"FBX6";"DKFZP564A2416";"FLJ10335";"STAT1";"FABP5";"MTHFD1";"CASP3" 

#cyan
moduleGenes[(GS > 0.4 & MM > 0.85)]
#"IFNG" "FIP2"


###【导出基因所在的模块】
probes = colnames(exp_mt_f)#【读取基因名】
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)#【读取模块】
geneInfo = geneInfo0[geneOrder, ]
#【保存所有模块基因】
write.table(geneInfo, file = "all_genes1.csv",sep="\t",row.names=F,quote=F)
###【输出每个模块的基因】
#【使用循环对每个模块的基因进行导出】
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(exp_mt_f)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}




#挑选模块Hub基因
## (1) 模块内基因连接度
adjacency = adjacency(exp_mt_f, power = 14)
TOM = TOMsimilarity(adjacency)
TOM[1:4,1:4]
Alldegrees =intramodularConnectivity(adjacency, net$colors)
head(Alldegrees)
#### kTotal:基因在整个网络中的连接度
#### kWithin: 基因在所属模块中的连接度，即Intramodular connectivity
#### kOut: kTotal-kWithin
#### kDiff: kIn-kOut

TOM <- as.matrix(TOM)
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
mycolor <- gplots::colorpanel(250, 'red', "orange", 'lemonchiffon')
TOMplot(
  plotTOM, net$dendrograms[[1]],   # 取对应块基因的颜色
  moduleColors[net$blockGenes[[1]]],
  main = "Network heatmap plot, all genes",
  col = mycolor
)

##绘制TOM热图
# tom plot
nGenes = ncol(exp_mt_f)
nSamples = nrow(exp_mt_f)
geneTree = net$dendrograms[[1]]
dissTOM = 1 - TOMsimilarityFromExpr(exp_mt_f, power = 6)
nSelect = 400
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected
genes",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))

##
MEs0 = moduleEigengenes(exp_mt_f, moduleColors)$eigengenes
# 计算根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs0);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
plotEigengeneNetworks(MEs0, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)

####----数据导出用于Cytoscape----####
#【构建文件夹】
cytoDir="CytoscapeInput"
dir.create(cytoDir)
for (mod in 1:nrow(table(moduleColors))){ #【对每个模块进行循环】
  modules = names(table(moduleColors))[mod]  
  probes = colnames(exp_mt_f)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes  
  modTOM = TOM[inModule, inModule] 
  dimnames(modTOM) = list(modProbes, modProbes)
  edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
  nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
  outEdge=paste(cytoDir,edges_File,sep="\\")
  outNode=paste(cytoDir,nodes_File,sep="\\")
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}
