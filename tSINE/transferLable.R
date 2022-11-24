rm(list=ls())
.libPaths(c('/local/txm/R/x86_64-pc-linux-gnu-library/4.0/',
            '/local/txm/R/x86_64-pc-linux-gnu-library/4.1/'))
setwd('/mdshare/node9/yanzijun/CRU/')
library(Seurat)
library(ggplot2)
library(scales)
library(reshape2)

########
## 导入单细胞数据
########
pbmc.Tcell <- readRDS('scTALL/public/NC_2021/single_cell_rnaseq_input/t_cell_atlas_subset.Rds')
# 确保RunUMAP加了return.model = T，便于最后整合在一个uMAP显示
pbmc_ref <- pbmc.Tcell
pbmc_ref <- FindVariableFeatures(object = pbmc_ref, selection.method = "vst", nfeatures = 3000)
pbmc_ref <- ScaleData(object = pbmc_ref, features = VariableFeatures(object = pbmc_ref))
pbmc_ref <- RunPCA(pbmc_ref, seed.use=123, npcs=50,
                   features = VariableFeatures(object = pbmc_ref), ndims.print=1,nfeatures.print=1)
pbmc_ref <- RunUMAP(pbmc_ref, dims = 1:50, seed.use = 123,n.components=2,
                      return.model = T)
Idents(pbmc_ref) <- pbmc_ref$Anno_level_2
DimPlot(pbmc_ref)

#######
### 导入查询数据
########
CL <- read.csv('TALL/data/SY/TALLcelline_RNAseq_Count.csv')
rownames(CL) <- CL$Symbol;CL$Symbol=NULL

## 构造seurat object
query=CL
pbmc_query <- CreateSeuratObject(counts= query, project = "cellline")
pbmc_query <- NormalizeData(object = pbmc_query, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc_query <- FindVariableFeatures(object = query, selection.method = "vst", nfeatures = 3000)
pbmc_query <- ScaleData(object = pbmc_query, features = VariableFeatures(object = pbmc_query))
pbmc_query <- RunPCA(object = pbmc_query, seed.use=123, npcs=50,
                     features = VariableFeatures(object = pbmc_query), ndims.print=1,nfeatures.print=1)
pbmc_query <- RunUMAP(pbmc_query, dims = 1:50, seed.use = 123,n.components=2,
                      return.model = T)
DimPlot(pbmc_query,reduction  = 'pca',label = T)
DimPlot(pbmc_query,label=T)

#######
### 转移标签
########
# 识别参考数据集的anchors
anchors <- FindTransferAnchors(reference = pbmc_ref, query = pbmc_query, dims = 1:30,
                               reference.reduction = 'pca')

# 将查询数据集映射到参考数据集上
predictions <- TransferData(anchorset = anchors,
                            refdata = pbmc_ref$Anno_level_2, dims = 1:30)
# 添加预测出的信息
pbmc_query <- AddMetaData(pbmc_query, metadata = predictions)
DimPlot(pbmc_query,group.by = 'predicted.id',reduction = 'pca')

pbmc_query$prediction.match <- pbmc_query$predicted.id == pbmc_query$Anno_level_2
table(pbmc_query$prediction.match)

pbmc_query <- MapQuery(anchorset = anchors, reference = pbmc_ref, query = pbmc_query, 
              refdata = list(celltype = "Anno_level_2"),reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(pbmc_ref, reduction = "umap", group.by = "Anno_level_2", label = TRUE, label.size = 3, 
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pbmc_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, 
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

