gset <- getGEO("GSE26713",destdir = '/local/yanzijun/CRU/TALL/data/GSE26713/')
meta <- pData(phenoData(gset[[1]]))
write.table(meta,'/local/yanzijun/CRU/TALL/data/GSE26713/metaInfo.txt',
             quote = F,sep='\t',row.names = T,col.names = T)
