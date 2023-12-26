library(affy)
setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE77939/rawdata")

# 读取数据
data <- ReadAffy(celfile.path = "C:/Users/11146/Desktop/双疾病生信数据/GSE77939/rawdata")

#把cdf文件转化为包，缺哪个注释文件就命名为哪个
#BiocManager::install('makecdfenv')
library(makecdfenv)
make.cdf.package("fatmito1a520158f.cdf", species = "Homo sapiens", compress = TRUE)
install.packages("fatmito1a520158fcdf", repos =  NULL, type="source" )
library(fatmito1a520158fcdf)

##归一化
eset.mas5 <- mas5(data)
eset.rma <- rma(data)

#提取表达矩阵
exprSet <- exprs(eset.mas5)
exprSet <- exprs(eset.rma)
exprSet <- log2(exprSet+1)
#更改列名
# 获取原始列名
colnames <- colnames(exprSet)
# 对所有列名进行替换
new_colnames <- gsub("\\.CEL\\.gz", "", colnames)
new_colnames <- sub("_.*", "", colnames)
# 将替换后的列名重新赋值给数据框
colnames(exprSet) <- new_colnames

save(exprSet,file = "exp.Rdata")
write.table(exprSet, "expr_mas5_matrix.txt", quote=F, sep="\t")
write.table(exprSet, "expr_rma_matrix.txt", quote=F, sep="\t")


