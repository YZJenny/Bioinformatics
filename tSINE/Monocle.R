library(monocle)
library('Seurat')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAST")path='/Users/yzj/Desktop/breast/data/29795293_Epith/'
ind4 <- read.table(paste(path,'GSE113196_RAW/GSM3099846_Ind4_Expression_Matrix.txt',sep=''),sep='\t', row.names=1, header=T)

pbmc <- CreateSeuratObject(counts = ind4, project = "IND", min.cells = 0, min.features = 0)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
dim(pbmc@assays$RNA)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#3.Feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#4.Scale
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#5.Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#6.Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc)
#7.tsne
pbmc_tsne <- RunTSNE(pbmc, dims = 1:20)
VEC <- pbmc_tsne@reductions$tsne@cell.embeddings
set.seed(12345)
KM=kmeans(VEC,centers=3)
KMC=KM$cluster
#plot(VEC,col=KMC,pch=15)
CLUST=KMC
CLUST=as.factor(CLUST)
names(CLUST)=names(pbmc_tsne@active.ident)
pbmc_tsne@active.ident=CLUST
DimPlot(pbmc_tsne, reduction = "tsne",label=T,pt.size=2)

pbmc.markers <- FindAllMarkers(pbmc_tsne, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE) #find markers for every cluster compared to all remaining cells, report only the positive ones
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) #plot the top 10 markers
DoHeatmap(pbmc_tsne, features = top10$gene) + NoLegend()
FeaturePlot(pbmc, cols=c("lightgrey",'red'), features = c("KRT18","KRT14","SLPI","ANKRD30A"))

####Monocle
count_matrix=as.matrix((pbmc_tsne@assays$RNA@counts))
pd=data.frame(sample=names(pbmc_tsne@active.ident),celltype=pbmc_tsne@active.ident)
pd1 <- new('AnnotatedDataFrame', data =pd)
fd=data.frame(gene_short_name=rownames(count_matrix))
rownames(fd)=rownames(count_matrix)
fd1 <- new('AnnotatedDataFrame', data =fd)

HSMM <- newCellDataSet(count_matrix,phenoData = pd1,featureData = fd1,expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

#Trajectory step 1: choose genes that define a cell's progress
ordering_genes <- as.character(read.table('/Users/yzj/Desktop/brain/Nguyen_Pervolarakis_Nat_Comm_2018/Monocle_Analysis/ind4_markergene.txt',header=T,sep='\t',row.names=1)[,6])

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

#Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "celltype")
