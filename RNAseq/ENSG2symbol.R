rm(list=ls())

## 处理gtf文件
awk '{if(!NF || /^#/){next}}1' Homo_sapiens.GRCh38.103.gtf | awk '{$2=null;$6=null;$7=null;$8=null;print $3"\t"$0}'| awk '{if(/^g/){print $0}}'| awk '{print $11"\t"$9"\t"$7"\t"$2"\t"$4"\t"$5}' |sed 's/"//g'| sed 's/;//g' > humanGTF

## ENSG2symbol
library(tidyverse)

count <- read.table('/remote-home/yanzijun/CRU/TALL/data/GSE227832/GSE227832_RNAseq_read_counts.txt',header = TRUE,sep='\t')
count[1:3,1:4]
colnames(count)[1] <- 'gene_id'

humanGTF <- read_delim("/remote-home/yanzijun/CRU/TALL/data/GSE227832/Homo_sapiens.GRCh38.103.txt", delim = "\t", 
                       escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE) %>% dplyr::select(X1,X3)
colnames(humanGTF) <- c("symbol","gene_id")
humanGTF$gene_id <- str_split(humanGTF$gene_id,"[.]",simplify = T)[,1] # 删除版本号，这里需要对ENSEMBL进行整理，删除“.”和后面的数字
humanGTF <- unique(humanGTF) %>% dplyr::select(gene_id,symbol) # 去重

humanGTF <- humanGTF[which(humanGTF$gene_id %in% count$gene_id),]
countf3 <- left_join(humanGTF,count,by="gene_id")
countf3 <- aggregate(x = countf3[,3:ncol(countf3)],   #此时exprSet的第三列开始是表达矩阵内容
                     by = list(symbol = countf3$symbol),   #按照相同symbol分组，在组内计算
                     FUN = mean) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')
write.table(countf3,'/remote-home/yanzijun/CRU/TALL/data/GSE227832/GSE227832_RNAseq_read_counts_symbol.txt',quote=F,row.names = F,col.names = TRUE)





data1 <- read.csv('/remote-home/yanzijun/CRU/TALL/data/GSE227832/GSE227832-GPL11154_series_matrix.txt',skip = 33,sep='\t')
data1$X.Sample_title <- NULL
data1 <- data1[,grep('diagnosis',colnames(data1))]
colnames(data1)


data2 <- read.csv('/remote-home/yanzijun/CRU/TALL/data/GSE227832/GSE227832-GPL16791_series_matrix.txt',skip = 33,sep='\t')
data2$X.Sample_title <- NULL
data1 <- data1[,grep('diagnosis',colnames(data1))]
colnames(data2)

data3 <- read.csv('/remote-home/yanzijun/CRU/TALL/data/GSE227832/GSE227832-GPL24676_series_matrix.txt',skip = 33,sep='\t')
data3$X.Sample_title <- NULL
data3 <- data3[,grep('diagnosis',colnames(data1))]
colnames(data2)
