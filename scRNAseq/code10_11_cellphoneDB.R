rm(list=ls())
setwd('/mdshare/node9/yanzijun/Extra/230307_code/')

library("stringr")
library(tibble)
library(Seurat)
library(ggplot2)
source('cellphone_util.R')
source('cellphone_dotplot.R')
source('cellphone_heat.R')
source('cellphone_cytoscape.R')
source('cellphone_circle.R')

pbmc <- readRDS('code4_pbmc.Rdata')

cells <- colnames(pbmc)[pbmc$type=='Microtia']
subpbmc <- subset(pbmc,cells = cells)

## 1. prepare data
data.input  <- as.data.frame(subpbmc@assays$RNA@data)
data.input <- rownames_to_column(data.input,var = "Gene")

meta = data.frame(Cell = colnames(subpbmc),CellType =subpbmc$subcelltype)
write.table(meta,'CellphoneDB_in/Microtia_meta.txt',row.names = FALSE,quote=FALSE,sep='\t')
write.table(data.input,'CellphoneDB_out/Microtia_count.txt',row.names = FALSE,quote=FALSE,sep='\t')

## 2. linux环境下跑 CellPhoneDB
# 安装
system('python -m venv cpdb-venv')
system('source cpdb-venv/bin/activate')
system('pip install cellphonedb')
# 运行
system('cellphonedb method statistical_analysis CellphoneDB_in/Microtia_meta.txt CellphoneDB_in/Microtia_count.txt --counts-data gene_name --threshold 0.1 --pvalue 0.05 --threads 20 --output-path CellphoneDB_out &')


## 3. 可视化
## code10 热图
heatmaps_plot(meta_file = 'CellphoneDB_out/Microtia_meta.txt',
                pvalues_file = 'pvalues.txt',count_filename = 'code10_heatmap_count.pdf',
                log_filename = 'code10_heatmap_log.pdf',
                show_rownames = T,show_colnames = T,scale="none",cluster_cols = F,border_color='white',cluster_rows = F,fontsize_row=20,
                fontsize_col = 20, main = '',treeheight_row=0,family='Arial',treeheight_col = 0,col1 = "dodgerblue4",col2 = 'peachpuff',col3 = 'deeppink4',
                meta_sep='\t',pvalues_sep='\t',pvalue=0.05)

# code11 点图
dot_plot(selected_rows = NULL,selected_columns = NULL,filename = 'code11_dotplot1.pdf',
         pvalue = 0.05,mean = 0,width = 30,height = 75,
         means_path = './CellphoneDB_out/means.txt',pvalues_path = './CellphoneDB_out/pvalues.txt', output_extension = '.pdf')
