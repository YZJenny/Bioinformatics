### bulk RNAseq数据的差异表达分析，输入数据类型包括COUNT/logTPM
rm(list=ls())
library(limma)
library(dplyr)
library(edgeR)
library(DESeq2)
library("BiocParallel") #启用多核计算

setwd('/Users/zijunyan/Desktop/临床中心/科研项目/T-ALL/data/2303_RNAseq_CQ/')

type <- 'COUNT' # logTPM,COUNT
TALL <- read.csv(paste('res/onlyTALL_Batch_after_gene_',type,'.csv',sep=''),row.names = 1)
load(paste('res/noTALL_Batch_after_gene_',type,'.RData',sep=''))

get_limma <- function(type,case.exp,ctrl.exp,case.tag,ctrl.tag,FC=1){
  ## 有两种方法:limma-trend和voom,在样本测序深度相差不大时两种方法差距不大，而测序深度相差大时voom更有优势，因此一般选择voom方法
  ## 准备数据
  sub.exp <- cbind(case.exp,ctrl.exp)
  group <- factor(c(rep(case.tag,ncol(case.exp)),rep(ctrl.tag,ncol(ctrl.exp)))) # levels顺序没关系

  ## 分组矩阵design构建
  design <- model.matrix(~0+group)
  rownames(design) = colnames(sub.exp)
  colnames(design) <- levels(group)
  
  ## 构建比较矩阵,比对顺序为case在前,control在后
  cont.matrix <- makeContrasts(contrasts = paste(case.tag,'-',ctrl.tag,sep=''), 
                               levels = design) ##比对顺序实验/对照!!!
  
  if(type=='COUNT'){
    ## 表达矩阵DGEList构建与过滤低表达基因
    dge <- DGEList(counts = sub.exp, group = group)
    keep.exprs <- filterByExpr(dge,design=design) #过滤低表达基因
    dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
    dge <- calcNormFactors(dge) #归一化基因表达分布
    
    ## DE分析:limma-trend(logCPM,有相似文库大小) or voom(文库大小差异大)
    # de <- cpm(dge, log=TRUE, prior.count=3)  #如选择logCPM，则eBayes设trend=TRUE
    de <- voom(dge, design, plot = FALSE, normalize = "quantile")
  }else{
    de  <- sub.exp
  }
  
  fit <- lmFit(de, design)
  fit1 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit1,trend=F)
  limma = topTable(fit2, coef = paste(case.tag,'-',ctrl.tag,sep=''), n = Inf) ##比对顺序实验/对照!!!
  limma = na.omit(limma)
  
  DEG <-limma %>% dplyr::filter((logFC> log(FC,2) | logFC < -log(FC,2) ) & adj.P.Val < 0.05)
  
  if(type=='COUNT'){
    write.csv(limma,paste('res/limma_COUNT/',case.tag,'vs',ctrl.tag,'.csv',sep=''),row.names = TRUE)
    write.csv(DEG,paste('res/limma_COUNT/',case.tag,'vs',ctrl.tag,'_FC',FC,'.csv',sep=''),row.names = TRUE)
  }else{
    write.csv(limma,paste('res/limma_logTPM/',case.tag,'vs',ctrl.tag,'.csv',sep=''),row.names = TRUE)
    write.csv(DEG,paste('res/limma_logTPM/',case.tag,'vs',ctrl.tag,'_FC',FC,'.csv',sep=''),row.names = TRUE)
  }
  
  res <- list(res=limma,DEG=DEG)
  return(limma)
}

get_edgeR <- function(case.exp,ctrl.exp,case.tag,ctrl.tag,FC=1){
  ## 准备数据,注意group需要relevel指定control
  sub.exp <- cbind(case.exp,ctrl.exp)
  group <- factor(c(rep(case.tag,ncol(case.exp)),rep(ctrl.tag,ncol(ctrl.exp))))
  group <- relevel(group,ctrl.tag) #将control组的因子设置为1!!!

  ## 分组矩阵design构建
  design <- model.matrix(~0+group)
  rownames(design) = colnames(sub.exp)
  colnames(design) <- levels(group)
  
  ## 表达矩阵DGEList构建与过滤低表达基因
  dge <- DGEList(counts = sub.exp,group = group) 
  keep.exprs <- filterByExpr(dge,design=design) #过滤低表达基因
  dge <- dge[keep.exprs,,keep.lib.sizes=FALSE] 
  dge <- calcNormFactors(dge, method = 'TMM') #归一化因子用于 normalizes the library sizes
  dge <- estimateDisp(dge, design, robust=T) 
  
  ##DE分析:官方建议bulk RNA-seq选择quasi-likelihood(QL) F-test tests;scRNA-seq 或是没有重复样品的数据选用likelihood ratio test
  fit <- glmQLFit(dge, design, robust=T)  #拟合模型 
  fit1 <- glmQLFTest(fit, contrast=c(-1,1)) #统计检验,注意比对顺序:实验-1 /对照1
  
  edgR <- topTags(fit1,n=Inf) #sort by p value, n is the number of genes/tags to return
  edgR <- na.omit(as.data.frame(edgR))
  
  DEG <- edgR %>% dplyr::filter((logFC> log(FC,2) | logFC < -log(FC,2) ) & FDR < 0.05)
  
  write.csv(edgR,paste('res/edgeR/',case.tag,'vs',ctrl.tag,'.csv',sep=''),row.names = TRUE)
  write.csv(DEG,paste('res/edgeR/',case.tag,'vs',ctrl.tag,'_FC',FC,'.csv',sep=''),row.names = TRUE)
  
  res <- list(res=edgR,DEG=DEG)
  return(edgR)
}

get_DESeq2 <- function(case.exp,ctrl.exp,case.tag,ctrl.tag,FC=1){
  ## 准备数据
  sub.exp <- cbind(case.exp,ctrl.exp)
  group <- factor(c(rep(case.tag,ncol(case.exp)),rep(ctrl.tag,ncol(ctrl.exp)))) #这里levels顺序没什么关系
  colData <- data.frame(row.names = colnames(sub.exp),
                        condition=group)
  
  ## 构建dds DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData=sub.exp,colData=colData,design=~condition)
  dds$condition <- relevel(dds$condition, ref = ctrl.tag)   #指定 control group
  keep <- rowSums(counts(dds)) >= 1.5*ncol(sub.exp)  #Pre-filtering ，过滤低表达基因
  dds <- dds[keep,] 
  
  ## DE分析
  dds <- DESeq(dds)
  DEseq <- results(dds,contrast=c('condition',case.tag,ctrl.tag)) #case在前,control在后，FC=case/control!!
  DEseq <- na.omit(as.data.frame(DEseq))
  
  DEG <- DEseq %>% dplyr::filter((log2FoldChange> log(FC,2) | log2FoldChange < -log(FC,2) ) & padj < 0.05)
  
  write.csv(DEseq,paste('res/DESeq2/',case.tag,'vs',ctrl.tag,'.csv',sep=''),row.names = TRUE)
  write.csv(DEG,paste('res/DESeq2/',case.tag,'vs',ctrl.tag,'_FC',FC,'.csv',sep=''),row.names = TRUE)
  
  res <- list(res=DEseq,DEG=DEG)
  return(DEseq)
}


## limma
limma_COUNT_vsThy <- get_limma(type='COUNT',case.exp=TALL,ctrl.exp=Thy,case.tag='TALL',ctrl.tag='Thy',FC=1)

limma_TPM_vsThy <- get_limma(type='logTPM',case.exp=TALL,ctrl.exp=Thy,case.tag='TALL',ctrl.tag='Thy',FC=1)

## edgR
edgR_vsThy <- get_edgeR(case.exp=TALL,ctrl.exp=Thy,case.tag='TALL',ctrl.tag='Thy',FC=1)

## DESeq2
DESeq2_vsThy <- get_DESeq2(case.exp=TALL,ctrl.exp=Thy,case.tag='TALL',ctrl.tag='Thy',FC=1)

### 4种差异分析结果比较
allg <- intersect(rownames(limma_COUNT_vsThy),rownames(edgR_vsThy))#取交集
allg <- intersect(allg,rownames(DESeq2_vsThy))
allg <- intersect(allg,rownames(limma_TPM_vsThy))
ALL_DEG <- cbind(limma_COUNT_vsThy[allg,c(1,4,5)],
                 limma_TPM_vsThy[allg,c(1,4,5)],
                 edgR_vsThy[allg,c(1,4,5)],
                 DESeq2_vsThy[allg,c(2,5,6)]) 
colnames(ALL_DEG)
colnames(ALL_DEG) <- c('limma_log2FC','limma_pvalue','limma_padj',
                       'limma_TPM_log2FC','limma_TPM_pvalue','limma_TPM_padj',
                       'edgeR_log2FC','edgeR_pvalue','edgeR_padj',
                       'DEseq2_log2FC','DEseq2_pvalue','DEseq2_padj')
##查看FC的相关性/一致性
print(cor(ALL_DEG[,c(1,4,7,10)]))

#查看显著差异基因重叠性，绘制韦恩图
#BiocManager::install("RBGL") #安装依赖包
#install.packages("Vennerable", repos="http://R-Forge.R-project.org") #安装Vennerable包
library(Vennerable) 

log2FC_cutoff=log2(1.5);  p_cutoff=0.05   #筛选显著差异基因比较
if(T){#根据Padj筛选
  DEseq2_deg <- rownames(DESeq2_vsThy[with(DESeq2_vsThy,abs(log2FoldChange)>log2FC_cutoff & padj<p_cutoff),])
  edgeR_deg <- rownames(edgR_vsThy[with(edgR_vsThy,abs(logFC)>log2FC_cutoff & FDR<p_cutoff),])
  limma_deg <- rownames(limma_COUNT_vsThy[with(limma_COUNT_vsThy,abs(logFC)>log2FC_cutoff & adj.P.Val<p_cutoff),])
}

mylist <- list(DEseq2=DEseq2_deg, edgeR=edgeR_deg, limma=limma_deg)
str(mylist)
Vennplot <- Venn(mylist)
Vennplot1 <- Vennplot[,c('DEseq2','edgeR')]
Vennplot2 <- Vennplot[,c('DEseq2','limma')]
Vennplot3 <- Vennplot[,c('limma','edgeR')]

pdf(file = paste0('fig/3DEG_Vennplot_lg2FC',log2FC_cutoff,'.pdf'))
plot(Vennplot, doWeights = F)
plot(Vennplot, doWeights = T) #doWeights=T设置为按数量比例绘图
plot(Vennplot1, doWeights = T)
plot(Vennplot2, doWeights = T)
plot(Vennplot3, doWeights = T)
dev.off()
