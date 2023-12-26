#genepix的芯片
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE19617/rawdata")
#提前解压好
# 指定所以后缀为.gpr文件的路径，方便read.maimages读取
files <- list.files("C:/Users/11146/Desktop/双疾病生信数据/GSE19617/rawdata/",
                    pattern = ".gpr")

path <- paste0("C:/Users/11146/Desktop/双疾病生信数据/GSE19617/rawdata/",
               files)

###########################读取文件###################
library(limma)
###干脆读取的时候把gene symbol一起提取出来
RG <- read.maimages(files=path,source="genepix") #读取文件
RG <- backgroundCorrect(RG, method="normexp") #背景校正
MA <- normalizeBetweenArrays(RG) #归一化
exprSet <- as.data.frame(MA$A) #获得表达矩阵
#修改列名
colnames <- colnames(exprSet)
new_colnames <- sub(".*/", "", colnames)
colnames(exprSet) <- new_colnames

#提取geneID
gene <- as.data.frame(MA$genes)
exp <- data.frame(gene$ID,exprSet) 


max(exprSet) #查看下最大值
boxplot(exprSet)#查看下整体分布
write.table(exp,"exp.csv",quote=F, sep="\t")
save(exp,file = "exp.Rdata")
