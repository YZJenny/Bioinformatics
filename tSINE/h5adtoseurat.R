h5adtoseurat <- function(file){
  library(reticulate)
  library(Seurat)
  ######加载python模块######
  scanpy <- import("scanpy")
  pandas <- import("pandas")
  
  prefix=unlist(strsplit(file,'\\.h5ad'))[1]
  adata = scanpy$read(file) ###载入scanpy输出的h5ad文件
  
  #######导出基因名和样本信息################
  meta = adata$obs
  gene <- adata$var
  
  #############导出矩阵并转置，scanpy和Seurat的行列是反的#############
  adata2 = adata$T
  adata2 = adata2$X
  
  ##给稀疏矩阵加个行名列名
  adata2@Dimnames[[1]] = rownames(gene)
  adata2@Dimnames[[2]] = rownames(meta)
  
  adata2 <- as(adata2, "CsparseMatrix")
  
  merge <- CreateSeuratObject(adata2)
  merge <- AddMetaData(merge, meta)
  saveRDS(merge,paste(prefix,'.RDS',sep=''))
}
