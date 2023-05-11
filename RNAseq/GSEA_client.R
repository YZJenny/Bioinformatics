rm(list=ls())
setwd('~/Desktop/Extra/230503_GSE136825/')

### 1. 制作表达谱输入文件，前两列为NAME，DESCRIPTION，后面接着样本
data <- read.csv('中间文件/exp.csv',row.names = 1)
data[1:3,1:3]
tmp <- data.frame(NAME=rownames(data),DESCRIPTION=rep('no',nrow(data)))
data <- cbind(tmp,data)
data[1:3,1:3]
write.table(data,'GSEA_input/expression.txt',quote = F,sep='\t',row.names = F,col.names = T)

### 2.制作表型文件,前两行手动添加
group <- readRDS('中间文件/group.rds')
group <- t(as.data.frame(group))
write.table(group,'GSEA_input/phenotype.txt',quote = F,sep='\t',row.names = F,col.names = F)

### 3. 挑选JAK_STAT通路
#grep 'JAK_STAT' KEGG.gmt > JAK_STAT.gmt

### 4.使用桌面版GSEA运行
