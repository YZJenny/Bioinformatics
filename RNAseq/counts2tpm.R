rm(list=ls())
setwd('/local/yanzijun/SY/data/221111/')
load('T-ALL-5batch.RData')

library(biomaRt)
mart=useMart('ensembl')
listDatasets(mart)
bmart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                          dataset = 'hsapiens_gene_ensembl',
                          host='www.ensembl.org') #grch38.p13版本
feature_ids <- rownames(adjusted_data)
length(feature_ids)
#listAttributes(bmart)
attributes=c(
  #'entrezgene_id',
  'hgnc_symbol','chromosome_name','start_position','end_position')
filters='hgnc_symbol'
#feature_info <- biomaRt::getBM(attributes = attributes,mart = bmart) 
#length(intersect(unique(feature_info$hgnc_symbol),rownames(adjusted_data))) #20332

feature_info <- biomaRt::getBM(attributes = attributes,filters = filters,values = feature_ids,mart = bmart)
dim(feature_info)

## 去除重复，只要染色体位置确定的
feature_info_full <- feature_info[!grepl('^CHR|^GL',feature_info$chromosome_name),]
dim(feature_info_full)

print(length(unique(feature_info_full$hgnc_symbol))) #其实还是有几十个重复的，忽略
feature_info_full <- feature_info_full[!duplicated(feature_info_full$hgnc_symbol),]
dim(feature_info_full)

mm <- match(feature_info_full$hgnc_symbol,feature_ids)
count <- adjusted_data[mm,]
dim(count)
print(all(rownames(count)==feature_info_full$hgnc_symbol))


## 计算gene有效长度
eff_length <- abs(feature_info_full$end_position-feature_info_full$start_position)
feature_info_full <- cbind(feature_info_full,eff_length)
eff_length2 <- feature_info_full[,c(1,5)]

x <- as.data.frame(count/eff_length2$eff_length)
tpm=t(t(x)/colSums(x))*1e6
logTPM <- log(tpm+1,2)
write.csv(logTPM,'logTPM_afterbatch.csv')
save.image('T-ALL-5batch.RData')
