#提取目的基因表达矩阵
setwd("C:/Users/11146/Desktop")
target_gene <- read.table("WGCNA+DEGs.txt")

setwd("C:/Users/11146/Desktop/双疾病生信数据/GSE33463")
load("exp.Rdata")
# 提取目的基因的基因表达矩阵
target_exp <- exp[target_gene$V1, ]
target_exp <- t(target_exp)

library(glmnet)
#构建lasso模型
x <- as.matrix(target_exp)
y <- group_list
fit <- glmnet(x, y, family = "binomial", alpha=1)

pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

#10倍交叉验证
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#筛选的特征基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
head(lassoGene)
write.table(lassoGene, file="LASSO-gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

#特征基因矩阵
rt <- t(target_exp)
lassoexp <- rt[lassoGene,,drop=F]
lassoexp <- as.data.frame(lassoexp)
write.table(lassoexp, file="LASSO-geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)



#SVM-REF
library(e1071)
library(kernlab)
library(caret)


Profile <- rfe(x=target_exp,
               y=as.numeric(as.factor(group_list)),
               sizes = c(2,4,6,8, seq(10,40,by=3)),
               rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
               methods="svmRadial")

plot(Profile,type=c("g", "o"))

SVMfeatureGenes=Profile$optVariables
head(SVMfeatureGenes)
write.table(SVMfeatureGenes, file="SVM-gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
rt1=t(target_exp)
svmexp=rt1[SVMfeatureGenes,,drop=F]
svmexp=as.data.frame(svmexp)
write.table(svmexp, file="SVM-geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)
