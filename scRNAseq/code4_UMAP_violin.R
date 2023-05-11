rm(list=ls())
setwd('/mdshare/node9/yanzijun/Extra/230307_code/')
library(Seurat)
library(dplyr)
library(ggplot2)
###################
## Seurat preprocessing: Normalization, scale and PCA,UMAP
###################
# pbmc <- readRDS('/mdshare/node9/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
# Idents(pbmc) <- pbmc$celltype
# demo <- subset(pbmc,downsample=200)
# saveRDS(demo,'code4_pbmc.Rdata')
# pbmc_batch <- readRDS('/mdshare/node9/yzj/JingMA_NEW/res/QC/ALL/RDS/PBMC_QC.RDS')
# demo2 <- subset(pbmc_batch,downsample=100)
# saveRDS(demo2,'code4_pbmc_QC.Rdata')

pbmc <- readRDS('code4_pbmc_QC.Rdata')
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 4000)
pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(object = pbmc, seed.use=123, npcs=150,
                     features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
pbmc <- RunTSNE(pbmc, dims = 1:50, seed.use = 123,n.components=2)
pbmc <- RunUMAP(pbmc, dims = 1:50, seed.use = 123,n.components=2)
DimPlot(pbmc)

##################
# 1. UMAP plot
##################
pbmc <- readRDS('code4_pbmc.Rdata')

CT <- c('CSPC', 'Chond', 'SSPC', 'SC','Immune','PVC','Endo')
Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")
names(Color) <- CT

pdf("code4_UMAP_CellType.pdf",width = 6,height = 6)
DimPlot(pbmc, group.by='celltype',label=T,pt.size = 0.1,
             cols = Color)
dev.off()

pdf("code4_UMAP_batch.pdf",width = 6,height = 6)
DimPlot(pbmc, group.by='batch', label=F,pt.size = 0.2)+
  theme(axis.text = element_text(size=20),
        panel.background=element_rect(fill='transparent', color='black',size = 1.5),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()


##################
# 2. violion plot
##################
marker.genes <- c('CDH5','CLDN5','PDGFRB','ACTA2','PTPRC','HLA-DRA','COL1A1','LUM','VCAN','ACAN','COL9A2','CYTL1','ELN','COL2A1','EGR1','HES1')
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
  df.sub <- data.frame(expvalue = pbmc@assays$RNA@data[gene,],
                       gene = rep(gene, ncol(pbmc@assays$RNA@data)),
                       celltype = pbmc$celltype)
  df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = CT)
color.cell <- c("#A6CEE3" ,"#1F78B4","#CAB2D6","#33A02C","#E31A1C","#FF7F00" ,"#6A3D9A")
plot.vln <- 
  ggplot(data = df.plot, aes(x = gene, y = expvalue, color = celltype, fill = celltype)) + 
  geom_violin(trim = T, scale = 'width') + 
  scale_color_manual(labels = CT,values = color.cell) + 
  scale_fill_manual(labels = CT,values = color.cell) + 
  facet_grid( ~ celltype) + 
  theme_classic() + coord_flip() +
  stat_summary(fun= mean, geom = "point",shape = 23, size = 2, color = "black") + 
  labs(x = 'Gene', y = 'Expression Level') + 
  theme(axis.text.y = element_text(size = 8, color = "black", face = 'italic'), 
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8,color = "black"), 
        strip.text.x = element_text(size = 8, color = "black"), legend.position = 'none')
plot.vln
ggsave("code4_Vln.pdf",plot.vln,height = 8, width = 16, units = 'cm')

##################
# 3. Feature plot
##################
marker.genes <- c('CDH5','CLDN5','PDGFRB','ACTA2')
pdf("code4_feature.pdf",width = 6,height = 6)
FeaturePlot(pbmc,marker.genes,ncol = 2)
dev.off()
