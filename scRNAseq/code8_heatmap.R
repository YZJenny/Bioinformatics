rm(list=ls())
.libPaths(c("/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.1",
            "/local/txm/R/x86_64-pc-linux-gnu-library/4.1"))
library(Seurat)
library(ggplot2)
library(monocle)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/mdshare/node9/yanzijun/Extra/230307_code/')

## 1. 查看TF的表达热图
TFs <- read.table('/mdshare/node9/yzj/publicData/Human_TF/TF_names_v_1.01.txt')[,1]
print(length(TFs))

HSMM <- readRDS('code7_monocle_output.RDS')

cds_subset <- HSMM[intersect(TFs,rownames(HSMM)),]
CT <- pData(HSMM)$celltype

pseudotime_de <- differentialGeneTest(cds_subset,cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
BEAM_res=BEAM(cds_subset,branch_point = 1,cores = 1)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]

saveRDS(BEAM_res, file = 'code8_BEAM.RDS')

pseudotime_TFs <- row.names(subset(BEAM_res,qval<1e-4))
print(length(pseudotime_TFs))
tmp1=plot_genes_branched_heatmap(cds_subset[pseudotime_TFs,],
                                 branch_point = 1,
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)

pdf("code8_branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()

### 2. 查看点基因的表达趋势
candiGenes <- c('CDK6','TUBB','TYMS')
cds_subset <- HSMM[candiGenes,]
CT <- pData(HSMM)$celltype
pdf('code8_pseu_gene.pdf',width = 8,height = 3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype",ncol = 3)+
  scale_fill_manual(values = getPalette[1:length(unique(pData(HSMM)$celltype))])+
  scale_color_manual(values = getPalette[1:length(unique(pData(HSMM)$celltype))])
dev.off()
