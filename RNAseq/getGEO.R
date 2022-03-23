gset <- getGEO("GSE26713",destdir = '/local/yanzijun/CRU/TALL/data/GSE26713/')
meta <- pData(phenoData(gset[[1]]))
write.table(meta,'/local/yanzijun/CRU/TALL/data/GSE26713/metaInfo.txt',
             quote = F,sep='\t',row.names = T,col.names = T)

exprSet<- exprs(gset[[1]])
exprSet[1:2,1:4]
## 获取探针与基因的对应关系
gpl <- getGEO('GPL570')
GPL=Table(gpl)
dim(GPL)

## 有的"Gene Symbol" 里面，有两个基因SYMBOL，因此我们需要拆分开后再继续注释
library("plyr")
GPL <- GPL[GPL$`Gene Symbol` != '', ]
dim(GPL)
split_gene <- strsplit( as.character( GPL$`Gene Symbol` ), " /// ")
probe2gene <- mapply( cbind, GPL[, 1], split_gene )
ID2gene <- as.data.frame(do.call(rbind, probe2gene))
head(ID2gene);dim(ID2gene)

exprSet <- exprSet[ID2gene$V1, ]
dim(exprSet)
exprSet[1:5, 1:6]


## 有一些基因对应多个探针或者一个探针对应多个基因，所以我们这里只保留表达量最大的探针
ID2gene$max <- apply(exprSet, 1, max)
ID2gene <- ID2gene[order(ID2gene$V2,
                         ID2gene$max,
                         decreasing = T), ]
ID2gene <- ID2gene[!duplicated(ID2gene$V2), ]
ID2gene <- ID2gene[!duplicated(ID2gene$V1), ]
head(ID2gene)
anno = ID2gene

## 最后只需要把行名调成一致即可
tmp = rownames(exprSet)%in% anno[,1]
exprSet = exprSet[tmp,]
dim(exprSet)

anno = anno[match(rownames(exprSet),anno$V1),]
dim(exprSet)
dim(anno)
rownames(exprSet) = anno$V2
dim(exprSet)
exprSet[1:5,1:5]
exprSet.uniq <- exprSet[!duplicated(rownames(exprSet)),]
write.table(exprSet.uniq,'CRU/TALL/data/GSE26713/exprSet.txt',quote = F,sep='\t',
            col.names = T,row.names = T)

