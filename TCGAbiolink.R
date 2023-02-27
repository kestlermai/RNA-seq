#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)#针对TCGAbiolinks包下载的数据进行整理

#查看可以下载的数据类型
getGDCprojects()$project_id

READ <- GDCquery(project = "TCGA-READ",data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "STAR - Counts")
GDCdownload(READ, method="api", directory = "./TCGA-READ-data")

expr <- GDCprepare(query=READ)#合并并创建一个对象储存数据
save("TCGA-READ-expr.Rdata")

#可更换其他数据格式，
exp <- expr@assays@data@listData$unstranded#原始count数
rownames(exp) <- expr@rowRanges$gene_id
colnames(exp) <- expr@colData@listData$barcode#TCGA的barcode就是TCGA-02-0001-01C-01D-0182-01
colnames(exp) <- substr(colnames(exp),start = 1, stop = 15)#取前十五位

exp <- as.data.frame(exp)
save(exp,file = "TCGA-READ-exp.Rdata")

#基因id注释
#https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files，下载gencode.v22.annotation.gtf.gz
#TCGA改版后用的是hg38
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)
setwd("C:/Users/11146/Desktop/R包/TCGA")

gene_info <- rtracklayer::import("gencode.v22.annotation.gtf.gz")#免去解压，直接读取
gene_info <- as.data.frame(gene_info)

#提取gene_id和gene_name中的两列的所有行做新的数据框
anno <- gene_info[,c("gene_id","gene_name")] %>% unique.data.frame()
rownames(anno) <- anno$gene_id

#intersect函数看两个数据集交集
combine_gene <- intersect(rownames(exp),rownames(anno))

exp <- exp[combine_gene,]
anno <- anno[combine_gene,]
exp <- data.frame(symbol=anno$gene_name, exp, check.names = F)
dim(exp)

lncRNA = c("3prime_overlapping_ncRNA", "antisense", 
             "bidirectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", 
             "non_coding", "processed_transcript", "sense_intronic" , "sense_overlapping")
table(gene_info$gene_type %in% lncRNA)

mRNA = c("protein_coding")
table(gene_info$gene_type %in% mRNA)

lnc_anno = gene_info[gene_info$gene_type %in% lncRNA,c("gene_name","gene_id","gene_type")]
mRNA_anno = gene_info[gene_info$gene_type == "protein_coding",c("gene_name","gene_id","gene_type")]
save(lnc_anno,mRNA_anno,file = "gtf.anno.Rdata")

#数据过滤，TCGA为了保证geneID在不同癌症的一致性，导致很多样本的表达量为0
exp_filtered <- exp[apply(exp, 1, function(x) sum(x > 1)> 89*0.5), ]#保留75%的样本基因数大于1；过滤的标准可以是保留50%，exp大于10的基因
dim(exp_filtered)

#分别取表达矩阵和mRNA/lncRNA_list的ENSID交集：
#mRNA交集：
mRNA <- intersect(rownames(exp_filtered),mRNA_anno$gene_stable_ID)

#lncRNA交集：
lncRNA <- intersect(rownames(exp_filtered),lncRNA_anno$gene_stable_ID)

#查看过滤后的mRNA和lncRNA数量(交集部分)：
length(mRNA)
length(lncRNA)

#mRNA表达矩阵：从表达矩阵中提取交集部分(mRNA)：
mRNA_exp <- exp_filtered[mRNA,]

#lncRNA表达矩阵：从表达矩阵中提取交集部分(lncRNA)：
lncRNA_exp <- exp_filtered[lncRNA,]

dim(mRNA_exp)
dim(lncRNA_exp)#行数和交集gene数相等，确认无误，表达矩阵拆分完成
#保存仅包含mRNA和lncRNA的两个独立表达矩阵
save(mRNA_exp,lncRNA_exp,mRNA_list,lncRNA_list,file = c('mRNA_lncRNA_ENSID_exp.Rdata'))



counts <- as.data.frame(assay(expr))#默认提取counts数据
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand"))#提取TPM数据
FPKM <- as.data.frame(assay(expr,i = "fpkm_unstrand"))#提取TPM数据
#Counts = "unstranded"；tpm = "tpm_unstrand"；fpkm = " fpkm_unstrand"

data=as.data.frame(rowRanges(expr))#获取其它信息数据
#这里面就包括注释以及编码、非编码等等信息

#信息组合
expr_count = cbind(gene_type=data$gene_type,gene_name=data$gene_name,counts)
expr_tpm = cbind(gene_type=data$gene_type,gene_name=data$gene_name,TPM)
expr_fpkm = cbind(gene_type=data$gene_type,gene_name=data$gene_name,FPKM)

write.csv(expr_count,'TCGA-READ_count.csv', quote = F)
write.csv(expr_tpm,'TCGA-READ_tpm.csv', quote = F)
write.csv(expr_fpkm,'TCGA-READ_fpkm.csv', quote = F)

# 1.临床信息获取
library(rjson)
library(tidyverse)
json <- jsonlite::fromJSON("metadata.cart.2022-11-08.json")
View(json)
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))
clinical <- read.delim('clinical.tsv',header = T)
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])
clinical_matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
clinical_matrix <- clinical_matrix[,-1]

# 2.临床信息整理
clinical_matrix$OS.time<-ifelse(clinical_matrix$vital_status=='Alive',
                                clinical_matrix$days_to_last_follow_up,
                                clinical_matrix$days_to_death)
#对于dead样本，overall survival采用day_to_death，对于alive样本，overall survival采用day_to_last_follow_up
clinical_matrix$OS.stat<-ifelse(clinical_matrix$vital_status=='Alive',0,1) #存活为0，死亡为1
cli<-data.frame(clinical_matrix$entity_submitter_id,clinical_matrix$age_at_index,clinical_matrix$gender,
                clinical_matrix$ajcc_pathologic_stage,clinical_matrix$OS.stat,clinical_matrix$OS.time)
colnames(cli) <- c("Id","age","gender","pathologic_stage","OS.stat","OS.time") 

write.csv(cli,'clinical.csv',quote = F)
