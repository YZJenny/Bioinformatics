rm(list=ls())
library(biomaRt)
setwd('/mdshare/node9/yanzijun/')

##################
### 1. COUNT转TPM
##################
COUNT <- read.table('/mdshare/node9/yanzijun/SY/data/GSE141140_all.counts.exp.txt',header = T,sep='\t')
rownames(COUNT) <- COUNT$Gene;COUNT$Gene=NULL

mart=useMart('ensembl')
listDatasets(mart)
bmart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                          dataset = 'hsapiens_gene_ensembl',
                          host='www.ensembl.org') #grch38.p13版本
feature_ids <- rownames(COUNT)
length(feature_ids)
#listAttributes(bmart)
attributes=c(
  #'entrezgene_id',
  'hgnc_symbol','chromosome_name','start_position','end_position')
filters='hgnc_symbol'
#feature_info <- biomaRt::getBM(attributes = attributes,mart = bmart) 
#length(intersect(unique(feature_info$hgnc_symbol),rownames(COUNT))) #20332

feature_info <- biomaRt::getBM(attributes = attributes,filters = filters,values = feature_ids,mart = bmart)
dim(feature_info)

## 去除重复，只要染色体位置确定的
feature_info_full <- feature_info[!grepl('^CHR|^GL',feature_info$chromosome_name),]
dim(feature_info_full)

print(length(unique(feature_info_full$hgnc_symbol)))
feature_info_full <- feature_info_full[!duplicated(feature_info_full$hgnc_symbol),]
dim(feature_info_full)

ol <- intersect(feature_info_full$hgnc_symbol,feature_ids)
feature_info_full <- feature_info_full[feature_info_full$hgnc_symbol %in% ol,]
dim(feature_info_full)
feature_info_full$start_position <- as.numeric(feature_info_full$start_position)
feature_info_full$end_position <- as.numeric(feature_info_full$end_position)

count <- as.data.frame(t(dplyr::select(as.data.frame(t(COUNT)),feature_info_full$hgnc_symbol)))
print(all(rownames(count)==feature_info_full$hgnc_symbol))


## 计算gene有效长度
eff_length <- abs(feature_info_full$end_position-feature_info_full$start_position)
feature_info_full <- cbind(feature_info_full,eff_length)
eff_length2 <- feature_info_full[,c(1,5)]

x <- as.data.frame(count/eff_length2$eff_length)
tpm=t(t(x)/colSums(x))*1e6
logTPM <- log(tpm+1,2)
write.table(logTPM,file="/mdshare/node9/yanzijun/SY/data/GSE141140.ltpm.txt",sep="\t",quote=F)
print(max(logTPM))
