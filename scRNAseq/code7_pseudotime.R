rm(list=ls())
.libPaths(c("/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.1",
            "/local/txm/R/x86_64-pc-linux-gnu-library/4.1"))

# tall <- readRDS('/mdshare/node9/yanzijun/CRU/scTALL/res/Seurat/RDS/pbmc_merge.RDS')
# tall <- subset(tall,downsample = 100)
#saveRDS(tall,'/mdshare/node9/yanzijun/Extra/230307_code/code7_monocle_input.RDS')

########
###  psedotime
########
setwd('/mdshare/node9/yanzijun/Extra/230307_code/')
library(Seurat)
library(ggplot2)
library(monocle)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

tall <- readRDS('code7_monocle_input.RDS')
data <- as(tall@assays$RNA@counts, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = tall@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

## 归一化
monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)

HSMM <- monocle
HSMM <- detectGenes(HSMM,min_expr =3 )
expressed_genes <- row.names(fData(HSMM))[fData(HSMM)$num_cells_expressed >= 30]
print(length(expressed_genes))
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~celltype",cores = 1)
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval,decreasing = FALSE)][1:3000]

HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2, reduction_method = "DDRTree")
HSMM <- orderCells(HSMM)

saveRDS(HSMM,'code7_monocle_output.RDS')

## 画图
pdf('monocle_merge.pdf',width = 6,height = 5)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",state_number_size=2,cell_size=0.8,cell_link_size=2)+
  theme(axis.line.x  = element_line(size=0.6, colour = "black"),
        axis.line.y  = element_line(size=0.6, colour = "black"),
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"))

plot_cell_trajectory(HSMM, color_by = "celltype",state_number_size=0.2,cell_size=0.8,
                     cell_link_size=2)+
  theme(axis.line.x  = element_line(size=0.6, colour = "black"),
        axis.line.y  = element_line(size=0.6, colour = "black"),
        axis.text = element_text(size = 15,colour = 'black'),
        axis.title = element_text(size = 15,colour = "black"))+
  scale_color_manual(values = getPalette[1:7])

plot_cell_trajectory(HSMM, color_by = "celltype",cell_size=0.8) +
  facet_wrap(~celltype, nrow = 2)+
  scale_color_manual(values = getPalette[1:7])+
  theme(axis.line.x  = element_line(size=0.4, colour = "black"),
        axis.line.y  = element_line(size=0.4, colour = "black"),
        axis.text = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 10,colour = "black"))

plot_cell_trajectory(HSMM, color_by = "State",state_number_size=2,cell_size=0.8,cell_link_size=2)+
  theme(axis.line.x  = element_line(size=0.6, colour = "black"),
        axis.line.y  = element_line(size=0.6, colour = "black"),
        axis.text = element_text(size = 15,colour = 'black'),
        axis.title = element_text(size = 15,colour = "black"))+
  scale_color_manual(values = getPalette[1:11])
dev.off()


