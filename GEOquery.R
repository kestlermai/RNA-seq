setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE30310")
#数据下载
#BiocManager::install("GEOquery")
library(GEOquery)
library(stringr)

eSetGSE30310 <- getGEO("GSE30310", 
                destdir = '.', 
                getGPL = F)

#(1)提取表达矩阵exp
exp1 <- exprs(eSetGSE30310[[1]])
#(2)提取临床信息
pd1 <- pData(eSetGSE30310[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd1),colnames(exp1));p
if(!p) exp1 = exp1[,match(rownames(pd1),colnames(exp1))]
#(4)提取芯片平台编号
gpl <- eSetGSE30310[[1]]@annotation

##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd1$`hiv status:ch1`,"positive"),"HIV","control")
group_list = factor(group_list,
                    levels = c("control","HIV"))

##提取gene symbol
#获得探针和基因名
GPL9392 <-getGEO('GPL9392',destdir =".")
GPL9392_anno <- Table(GPL9392)
colnames(GPL9392_anno)

#提取对应gene symbol
anno <- GPL9392_anno[,c("ID","GeneName")]

load("exp.Rdata")
save(exprSet,group_list,file = "exp.Rdata")

library(tidyverse)
library(dplyr)
library(tibble)
exprSet <- as.data.frame(exp1)
exprSet <- exprSet %>% 
  mutate(ID=rownames(exprSet)) %>% 
  inner_join(anno,by="ID") %>%
  select(-ID) %>% #去掉多余信息
  select(GeneName, everything()) %>% #重新排列，
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(GeneName != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(GeneName,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "GeneName")
#tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除，变成行名后genename好像会丢失
write.csv(exprSet,"exp_raw.csv",row.names = T,quote = F)
#对GSE30310进行归一处理
exprSet <- exprSet[,1:67]
pd1 <- pd1[1:67,]
p = identical(rownames(pd1),colnames(exprSet));p
if(!p) exp1 = exp1[,match(rownames(pd1),colnames(exprSet))]
group_list=ifelse(str_detect(pd1$`hiv status:ch1`,"positive"),"HIV","control")
group_list = factor(group_list,
                    levels = c("control","HIV"))

library(limma)
exp_norm <- normalizeBetweenArrays(exprSet)
exp_norm -> exp
save(exp,group_list,file = "exp.Rdata")
write.csv(exp_norm,"exp_norm.csv",row.names = T,quote = F)

#绘制质量箱式图
png("boxplot.png", width=800, height=500)
#par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE9392")
boxplot(exp_norm, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()


#用geneExpressionFromGEO注释GSE131793和GSE703
#BiocManager::install("geneExpressionFromGEO")
#GPL11532，GPL23126，GPL6244，GPL8300，GPL80，GPL96，GPL570，GPL571，GPL2115，GPL1293，GPL6102
#GPL6104，GP16883，GPL6884，GP113497，GPL14550，GP117077，GP16480。
library(GEOquery)
library(geneExpressionFromGEO)
library(limma)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(cowplot)

setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE703")

library(GEOquery)
library(stringr)
eSet <- getGEO("GSE703", 
               destdir = '.', 
               getGPL = F)


exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$title,"PAH"),"PAH","control")
group_list = factor(group_list,
                    levels = c("control","PAH"))
##提取gene symbol
#获得探针和基因名
GPL80 <-getGEO('GPL80',destdir =".")
GPL80_anno <- Table(GPL80)
colnames(GPL80_anno)

#提取对应gene symbol
anno <- GPL80_anno[,c("ID","Gene Symbol")]
colnames(anno)[colnames(anno) == "Gene Symbol"] <- "symbol"
# 提取列中的第一个基因名称
anno$symbol <- sapply(strsplit(anno$symbol, "///"), "[", 1)

library(tidyverse)
library(dplyr)
library(tibble)
exprSet <- as.data.frame(exp)
exprSet <- exprSet %>% 
  mutate(ID=rownames(exprSet)) %>% 
  inner_join(anno,by="ID") %>%
  select(-ID) %>% #去掉多余信息
  select(symbol, everything()) %>% #重新排列，
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  filter(symbol != "---") %>%
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exprSet
exp <- log2(exp+1)
save(exp,group_list,file = "exp.Rdata")


#绘制质量箱式图
pdf("boxplot.pdf", width=8, height=5)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE703")
boxplot(exp, boxwex=0.9,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        col=c("#0072B5FF","#BC3C29FF"),
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()

#标化
exp <- normalizeBetweenArrays(exp)
save(exp,group_list,file = "exp.Rdata")
#ggplot2绘制
df1 <- melt(uniExpr)

# 密度图
densi <- ggplot(df1,aes(x = log2(value),color = variable)) +
  geom_density(linewidth = 1,show.legend = F) +
  theme_bw(base_size = 18)
densi
# 箱线图
box <- ggplot(df1,aes(x = variable,y = log2(value))) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('')
box
# 拼图
plot_grid(plotlist = list(densi,box),ncol = 2,align = 'hv')
ggsave("densi-box.png",width = 20, height = 10, dpi = 600)

#判断原始数据是否去了log
max(uniExpr)
if(max(uniExpr)>30) uniExpr=log2(uniExpr+1) #最大值大于30则取log

save(uniExpr,group_list,file = "exp.Rdata")

write.csv(uniExpr,"exp.csv",row.names = T,quote = F)

#GSE131793
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE131793")
library(GEOquery)
library(stringr)
eSet <- getGEO("GSE131793", 
               destdir = '.', 
               getGPL = F)

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$title,"Control"),"control","PAH")
group_list = factor(group_list,
                    levels = c("control","PAH"))
#GPL6244对应的是：hugene10sttranscriptcluster
BiocManager::install("hugene10sttranscriptcluster.db")
library("hugene10sttranscriptcluster.db")
ids <- toTable(hugene10sttranscriptclusterSYMBOL)

#去除空值symbol
ids<- ids[ids$symbol!=" ",]
ids <- ids[ids$probe_id %in%  rownames(exp),]
exp <- exp[rownames(exp) %in%ids$probe_id,  ]

#由于多个探针对应一个基因，所以取基因在样本中的均值并排序取最大得作为唯一探针与基因对应。
##加载注释函数
anno <- function(exp,ids){
  probes <- as.character(by(exp,ids$symbol,function(x)rownames(x)[which.max(rowMeans(x))]))
  exp=exp[rownames(exp) %in% probes ,]
  rownames(exp)=ids[ids$probe_id%in%rownames(exp),2]
  return(exp)
}
exp<- anno(exp,ids)
#保存注释后的表达矩阵以及分组信息
save(exp,group_list,file = "exp.Rdata")


setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE140713")
#数据下载
#BiocManager::install("GEOquery")
library(GEOquery)
library(stringr)

eSet<- getGEO("GSE140713", 
              destdir = '.', 
              getGPL = F)

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$`disase state:ch1`,"AIDS"),"HIV","control")
group_list = factor(group_list,
                    levels = c("control","HIV"))

library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL6480')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  mutate(probe_id=rownames(exp)) %>% 
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列

save(exp,group_list,file = "exp.Rdata")

pdf("boxplot.pdf", width=10, height=5)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE140713")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        col=c("#0072B5FF","#BC3C29FF"),
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()
library(limma)
exp <- normalizeBetweenArrays(exp)#用了会减少组内差异
save(exp,group_list,file = "exp.Rdata")


setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
#数据下载
#BiocManager::install("GEOquery")
library(GEOquery)
library(stringr)

eSet<- getGEO("GSE33463", 
              destdir = '.', 
              getGPL = F)

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$`disease status:ch1`,"control"),"control","PAH")
group_list = factor(group_list,
                    levels = c("control","PAH"))

library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL6947')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  mutate(probe_id=rownames(exp)) %>% 
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列

exp <- exp[,1:113]
pd <- pd[1:113,]
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
group_list=ifelse(str_detect(pd$`disease status:ch1`,"control"),"control","PAH")
group_list = factor(group_list,
                    levels = c("control","PAH"))
exp <- normalizeBetweenArrays(exp)
save(exp,group_list,file = "exp.Rdata")

pdf("boxplot.pdf", width=10, height=5)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE33463")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        col=c("#0072B5FF","#BC3C29FF"),
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()

setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE19617")
library(GEOquery)
library(stringr)
eSet <- getGEO("GSE19617", 
                destdir = '.', 
                getGPL = F)

exp <- read.delim("exp.txt",header = T,sep = "\t")
pd <- pData(eSet[[1]])
gpl <- eSet[[1]]@annotation

library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL6480')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列

#提取
exp1 <- exp[,-(1:25)]
pd1 <- pd[-(1:25),]
p1 = identical(rownames(pd1),colnames(exp1));p1
if(!p1) exp1 = exp1[,match(rownames(pd1),colnames(exp1))]
library(stringr)
group_list=ifelse(str_detect(pd1$title,"Normal"),"control","PAH")
group_list = factor(group_list,
                    levels = c("control","PAH"))
table(group_list)

exp <- normalizeBetweenArrays(exp)
save(exp,group_list,file = "exp.Rdata")


setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE77939")
library(GEOquery)
library(stringr)
eSet <- getGEO("GSE77939", 
               destdir = '.', 
               getGPL = F)
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE77939/rawdata")
load("exp.Rdata")
exp <- exprSet
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$`subject status:ch1`,"Healthy"),"control","HIV")
group_list = factor(group_list,
                    levels = c("control","HIV"))

##提取gene symbol
#获得探针和基因名
GPL15207 <-getGEO('GPL15207',destdir =".")
GPL15207_anno <- Table(GPL15207)
colnames(GPL15207_anno)

#提取对应gene symbol
anno <- GPL15207_anno[,c("ID","Gene Symbol")]
colnames(anno)[colnames(anno) == "Gene Symbol"] <- "symbol"
# 提取列中的第一个基因名称
anno$symbol <- sapply(strsplit(anno$symbol, "///"), "[", 1)

library(tidyverse)
library(dplyr)
library(tibble)
exprSet <- as.data.frame(exp)
exprSet <- exprSet %>% 
  mutate(ID=rownames(exprSet)) %>% 
  inner_join(anno,by="ID") %>%
  select(-ID) %>% #去掉多余信息
  select(symbol, everything()) %>% #重新排列，
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  filter(symbol != "---") %>%
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exprSet
exp <- normalizeBetweenArrays(exp)
save(exp,group_list,file = "exp.Rdata")

pdf("boxplot1.pdf", width=8, height=5)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE77939")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        col=c("#0072B5FF","#BC3C29FF"),
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()

setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33877")
library(GEOquery)
library(stringr)
eSet <- getGEO("GSE33877", 
               destdir = '.', 
               getGPL = F)

exp <- read.delim("norexp1.txt",header=T,sep = "\t")
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$`infection:ch1`,"HIV"),"HIV","control")
group_list = factor(group_list,
                    levels = c("control","HIV"))

#devtools::install_github("jmzeng1314/idmap3")
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL6947')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列
save(exp,group_list,file = "exp.Rdata")

png("boxplot.png", width=800, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE33877")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()



setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE27393")
#数据下载
#BiocManager::install("GEOquery")
library(GEOquery)
library(stringr)

eSet<- getGEO("GSE27393", 
               destdir = '.', 
               getGPL = F)

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
##创建分组信息
library(stringr)
group_list=ifelse(str_detect(pd$`strain:ch1`,"HIV"),"HIV","control")
group_list = factor(group_list,levels = c("control","HIV"))
#注释基因id
##提取gene symbol
#获得探针和基因名
GPL1355 <-getGEO('GPL1355',destdir =".")
GPL1355_anno <- Table(GPL1355)
colnames(GPL1355_anno)

#提取对应gene symbol
anno <- GPL1355_anno[,c("ID","Gene Symbol")]
# 提取列中的第一个基因名称
anno$`Gene Symbol` <- sapply(strsplit(anno$`Gene Symbol`, "///"), "[", 1)
#去除anno中空值
anno %>% filter(`Gene Symbol` != "")
  
library(tidyverse)
library(dplyr)
library(tibble)
exprSet <- as.data.frame(exp)
exprSet <- exprSet %>% 
  mutate(ID=rownames(exprSet)) %>% 
  inner_join(anno,by="ID") %>%
  select(-ID) %>% #去掉多余信息
  select(`Gene Symbol`, everything()) %>% #重新排列，
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(`Gene Symbol` != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(`Gene Symbol`,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "Gene Symbol")

exp <- log2(exprSet+1)
save(exp,group_list,file = "exp.Rdata")

png("boxplot.png", width=800, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE27393")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()


setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE68563")
eSet <- getGEO("GSE68563", 
               destdir = '.', 
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
#注释基因id
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL13497')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  mutate(probe_id=rownames(exp)) %>%
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列
save(exp,group_list,file = "exp.Rdata")

exp <- exp[,!colnames(exp) %in% c("GSM1675623","GSM1675626","GSM1675629")]

pd <- pd[!rownames(pd) %in% c("GSM1675623","GSM1675626","GSM1675629"),]

p = identical(rownames(pd),colnames(exp));p
##创建分组信息
library(stringr)
group_list <- ifelse(str_detect(pd$source_name_ch1,"HIV"),"HIV","control")
group_list <- factor(group_list,levels = c("control","HIV"))

exp <- log2(exprSet+1)
save(exp,group_list,file = "exp.Rdata")

png("boxplot.png", width=800, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE68563")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()


setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE58994")
eSet <- getGEO("GSE58994", 
               destdir = '.', 
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
#注释基因id
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL16686')
head(ids) 

GPL16686 <- getGEO('GPL16686',destdir =".")
GPL16686_anno <- Table(GPL16686)
colnames(GPL16686_anno)
library(org.Hs.eg.db)
library(clusterProfiler)
Gene_Anno <- bitr(GPL16686_anno$GB_ACC, fromType = "REFSEQ", #fromType是指你的数据ID类型是属于哪一类的
                  toType = "SYMBOL", #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = org.Hs.eg.db)#Orgdb是指对应的注释包是哪个
Gene_Anno
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  mutate(probe_id=rownames(exp)) %>%
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列
save(exp,group_list,file = "exp.Rdata")

colnames(exp)

exp <- exp[,!colnames(exp) %in% c("GSM1675623","GSM1675626","GSM1675629")]

pd <- pd[!rownames(pd) %in% c("GSM1675623","GSM1675626","GSM1675629"),]

p = identical(rownames(pd),colnames(exp));p
##创建分组信息
library(stringr)
group_list <- ifelse(str_detect(pd$source_name_ch1,"HIV"),"HIV","control")
group_list <- factor(group_list,levels = c("control","HIV"))

exp <- log2(exprSet+1)
save(exp,group_list,file = "exp.Rdata")

png("boxplot.png", width=800, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE68563")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()


setwd("C:/Users/maihuanzhuo/Desktop/双疾病生信数据/GSE33879")
eSet <- getGEO("GSE33879", 
               destdir = '.', 
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[2]])
#(2)提取临床信息
pd <- pData(eSet[[2]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[2]]@annotation
#注释基因id
library(idmap3)
ids <- idmap3::get_pipe_IDs('GPL6947')
head(ids) 
###id转换
library(tidyverse)
exp <- as.data.frame(exp) #也可以直接加到管道操作中
#id匹配
exp <- exp %>% 
  mutate(probe_id=rownames(exp)) %>%
  inner_join(ids,by="probe_id") %>%
  select(probe_id, symbol, everything()) %>%
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  filter(symbol != "") %>% #去除symbol中的空值
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  select(-rowMean) %>% #反向选择去除rowMean这一列
  column_to_rownames(var = "symbol")
exp <- exp[,-(1)] #去掉第一列
##创建分组信息
library(stringr)
group_list <- ifelse(str_detect(pd$`infection:ch1`,"HIV"),"HIV","control")
group_list <- factor(group_list,levels = c("control","HIV"))
exp <- log2(exp+1)
save(exp,group_list,file = "exp.Rdata")

png("boxplot.png", width=800, height=500)
par(mar=c(7,4,2,1))#底部、左侧、顶部和右侧的边距大小
title <- paste ("GSE33879")
boxplot(exp, boxwex=0.7,#箱子宽度
        notch=T,#notch表示是否在箱子上添加缺口
        main=title, cex.main=2,
        outline=FALSE, #outline确定是否显示异常值
        las=2)#las设置了轴标签的方向（2表示垂直）
dev.off()
