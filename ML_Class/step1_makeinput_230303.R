rm(list=ls())
library(dplyr)

setwd('/mdshare/node9/yanzijun/Extra/230228_Model/')
Expr <- read.table('data_mrna_seq_v2_rsem.txt',header = T,sep='\t')
Expr <- Expr[!duplicated(Expr$Hugo_Symbol),]
Expr <- na.omit(Expr)
rownames(Expr) <- Expr$Hugo_Symbol;
Expr$Hugo_Symbol=Expr$Entrez_Gene_Id=NULL
Expr[1:3,1:3]

## change colnames of data
change_colname <- function(old_col){
  new_col <- apply(as.matrix(old_col),1,function(x) 
    paste(unlist(strsplit(x,'\\.'))[1:3],collapse ='-')
  )
  return(new_col)
}

old_col <- colnames(Expr)
new_col <- change_colname(old_col)
colnames(Expr) <- new_col

#all.label <- read.csv('方法1_两组筛差异基因_N+ vs N0.csv')
#all.label <- read.csv('方法2_两组筛差异基因_N2-3 vs N0.csv')
all.label <- read.csv('方法3_两组筛差异基因_N3 vs N0.csv')
set.seed(123)
index <- sample(1:nrow(all.label),floor(nrow(all.label)*0.7))
train.label <- all.label[index,]
test.label <- all.label[-index,]
table(all.label$Group)
table(train.label$Group)
table(test.label$Group)

## 统一既有表达值又有临床数据的样本
train.label <- train.label[train.label$sample %in% intersect(colnames(Expr),train.label$sample),]
rownames(train.label) <- train.label$sample;train.label$sample=NULL

test.label <- test.label[test.label$sample %in% intersect(colnames(Expr),test.label$sample),]
rownames(test.label) <- test.label$sample;test.label$sample=NULL

train.exp <-select(Expr,rownames(train.label))
test.exp <-select(Expr,rownames(test.label))

print(all(colnames(train.exp)==rownames(train.label)))
print(all(colnames(test.exp)==rownames(test.label)))

write.table(train.exp,'input_method3/train_exp.txt',quote = F,sep='\t',col.names = T,row.names = T)
write.table(test.exp,'input_method3/test_exp.txt',quote = F,sep='\t',col.names = T,row.names = T)

write.table(train.label,'input_method3/train_label.txt',quote = F,sep='\t',col.names = T,row.names = T)
write.table(test.label,'input_method3/test_label.txt',quote = F,sep='\t',col.names = T,row.names = T)

data.lst <- list(train.exp=train.exp,test.exp=test.exp,train.label=train.label,test.label=test.label)
saveRDS(data.lst,'input_method3/data.lst')
