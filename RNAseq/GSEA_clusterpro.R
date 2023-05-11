rm(list = ls())
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

setwd('/mdshare/node9/yanzijun/Extra/230305_gsea/')

gene_FC <- read.csv('gsea_input.csv')
colnames(gene_FC) <- c('Gene','logFC')
gene_FC$Gene <- toupper(gene_FC$Gene)
genename <- gene_FC$Gene

## 1.将基因名转换成ENTREZID，如果用GSEABase包，就不需要转
# gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
# colnames(gene_map)[1]<-"Gene"
# write.csv(as.data.frame(gene_map),"基因ID转换.csv",row.names =F)#导出结果至默认路径下

## 2.将ENTREZID与logFC结合，并根据logFC的值降序排列
# aaa<-inner_join(gene_map,gene_FC,by = "Gene")
# aaa<-aaa[,-1]
# aaa<-na.omit(aaa)
# aaa$logFC<-sort(aaa$logFC,decreasing = T)

## 3.整理成GSEA分析的格式
geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList
#GSEA分析——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#保存文件
write.table (Go_gseresult, file ="Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="Go_Reactomeresult.csv", sep =",", row.names =TRUE)

## 5.画图
gseaplot2(egmt,1,pvalue_table = TRUE)#输出第212个结果
