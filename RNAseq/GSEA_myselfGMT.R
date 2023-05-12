rm(list = ls())
library(GSEABase) 
library(clusterProfiler)
setwd('/mdshare/node9/yanzijun/Extra/230305_gsea/')

gene_FC <- read.csv('gsea_input.csv')
colnames(gene_FC) <- c('Gene','logFC')
gene_FC$Gene <- toupper(gene_FC$Gene)

## 1.整理成GSEA分析的格式
gene_FC <- gene_FC[order(gene_FC$logFC,decreasing = T),]
geneList = gene_FC[,2]
names(geneList) = as.character(gene_FC[,1])
head(geneList);tail(geneList)

## 2.自己制作gmt文件
FG <- readRDS('/mdshare/node9/yanzijun/public/FerrDB/FC_dr_mk_sup.RDS')
gset <- c("ferroptosis_geneset","NA",FG)
gset <- gset%>% 
  as.data.frame() %>% 
  t()
write.table(gset,file = "ferroptosis_geneset.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

## 3. 分析
geneset <- read.gmt('ferroptosis_geneset.gmt')  
egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff = 1,seed = 123)
head(egmt)
write.table(egmt, file ="ferroptosis.csv", sep =",", row.names =FALSE)

## 4.画图
plot_gsea <- gseaplot2(egmt,1,pvalue_table = TRUE)
pdf('ferroptosis_gsea.pdf',plot_gsea,width = 7,height = 5)
print(plot_gsea)
dev.off()
