rm(list=ls())
setwd('/mdshare/node9/yanzijun/Extra/230307_code/')
library(Seurat)
library(dplyr)
library(ggplot2)

pbmc <- readRDS('code4_pbmc.Rdata')
Idents(pbmc) <- pbmc$celltype

CT <- c('CSPC', 'Chond', 'SSPC', 'SC','Immune','PVC','Endo')
Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")
names(Color) <- CT

### 1.2 all samples ratio in each cluster
phylog_df <- pbmc@meta.data[,c('batch',"celltype")]
phylog_df <- table(phylog_df$batch,phylog_df[,"celltype"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType,levels = rev(CT))

p <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.8)+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.title = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.title = element_text(size=15,face = 'bold',colour = 'black'),
        legend.text = element_text(size=15,face = 'bold',colour = 'black'))+
  labs(x='',y='Cell proportion')+theme(legend.position="right")+
  scale_fill_discrete(guide = guide_legend(reverse=TRUE))+
  scale_fill_manual(values = rev(Color))
p
ggsave('code5_cellProportion.pdf',p,width = 10,height = 10)