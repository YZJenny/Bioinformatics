rm(list=ls())
library(MuSiC)
library(Biobase)
library(Seurat)
setwd('/mdshare/node9/yanzijun/Extra/230225_MuSiC/')
bulk.eset = readRDS('s1407717_expression_set.rds')
bulk.mtx=bulk.eset@assayData$exprs
#原先bulk样本只有一例，会报错，所以copy一列，防止报错
bulk.mtx=as.data.frame(bulk.mtx)
bulk.mtx$s1407717_copy <- bulk.mtx$s1407717
bulk.mtx=as.matrix(bulk.mtx)

load('humanSCanno.RData') 
#seurat添加一列sampleID，并转成singlecellexperiment
humanP$sampleID=colnames(humanP)
sc.eset=as.SingleCellExperiment(humanP)

#remotes::install_github("renozao/xbioc")
library(xbioc)
library(SummarizedExperiment)
# Estimate cell type proportions

Est.prop = music_prop(bulk.mtx = bulk.mtx, sc.sce = sc.eset, 
                      markers = intersect(rownames(humanP),rownames(bulk.mtx)),
                               clusters = 'res0.2.names',samples = 'sampleID', 
                               verbose = F)
save.image('MuSiC.RData')

# plot estimated cell type proportions
print(names(Est.prop))
df <- as.data.frame(Est.prop$Est.prop.weighted[1,])
colnames(df)[1] <- colnames(bulk.mtx)[1]
df <- tibble::rownames_to_column(df,'celltype')

plot.MuSiC <- ggplot(df, aes(x = celltype,y=s1407717,fill=celltype))+
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = 'none')+ylab("Proportion")
ggsave('MuSiC.pdf',plot = plot.MuSiC,width = 5,height = 4)
