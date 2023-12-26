#illumina平台文件表达矩阵整理
library(limma)
setwd("C:/Users/11146/Desktop")

x <- read.ilmn("GSE33877_non-normalized.txt", 
               expr="AVG_Signal", #命名分割规律
               probeid = "ID_REF")
y <- neqc(x)
exp=y$E


rt2=rbind(ID=colnames(exp),exp)
write.table(rt2,file="norexp.txt",sep="\t",quote=F,col.names = F)


cols=rainbow(ncol(exp)) 
pdf(file = "nor.pdf",width=6,height = 4.5)
par(cex = 0.3,mar=c(8,8,8,8))
boxplot(exp,col=cols)
dev.off()

