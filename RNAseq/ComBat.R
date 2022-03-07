### ComBat消除batch,Combat的输入是标准化后的数据, 不能存在全0的行，以及Inf的出现
library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(sva)
library(dplyr)
library(edgeR)

zero.rows.del <- function(data,batch){
  data <- as.matrix(data)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(data[, batch == batch_level], 1,
                         function(x) {
                           var(x) == 0
                         })))
    }else {
      
      return(which(rep(1, 3) == 2))
      
    }       
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(data), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n",
                length(zero.rows)))
    data.orig <- data
    data <- data[keep.rows, ]
  }
  return(data)
}

  
counts1 <- read.table('/local/yanzijun/SY/data/GSE141140_all.counts.exp.txt',header = TRUE)
counts2 <- read.csv('/local/yanzijun/SY/data/AllSamples.GeneExpression.Count.csv',header = TRUE)[,c(1:55)]
counts2$gene_id <- NULL
colnames(counts2)[54] <- colnames(counts1)[1]
## counts1: 17 samples,counts2: 53 samples

counts <- merge(counts1,counts2,by='Gene')
rownames(counts) <- counts$Gene
counts$Gene <- NULL
dim(counts) #70

## Combat的输入是标准化后的数据，所以先进行cpm转换
library(edgeR)
lcpm <- cpm(counts, log=TRUE, prior.count=2)

## 构建表型信息
tissue <- rep('TALL',ncol(counts))
tissue[grep('^N',colnames(counts))] <- 'Tcell'
tissue[grep('^HSC',colnames(counts))] <- 'HSC'
tissue[grep('^Thy',colnames(counts))] <- 'Thy'
# tissue <- factor(tissue,levels = c('TALL','HSC','Thy','Tcell'))

batch <- paste0("batch", rep(c(1,2),c(17,53)))
table(batch,tissue) 

## 执行Combat
mod <- model.matrix(~tissue)
expr_batch <- ComBat(dat = lcpm, batch = batch, mod = mod)

## 检验效果
library("FactoMineR")
library("ggplot2")
library("factoextra")
pca.plot = function(dat,col){
  
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}

ind <- rep(c("Batch1","Batch2"),c(17,53))
pdf('SY/res/Combat.pdf',width = 5,height = 4)
print(pca.plot(lcpm,factor(ind)))  ## 处理前
print(pca.plot(expr_batch,factor(ind)))  ## 处理后
dev.off()


dat <- melt(lcpm)
p <- ggplot(dat, aes(x=Var2, y=value,color = Var2)) + 
  geom_boxplot() +
  theme_classic()+labs(title="Distribution of expression",y='exp',x='samples')+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5),
        legend.position="none")
ggsave('SY/res/Samples.pdf',p,width = 20,height = 4)


## 查看S1PR3的表达
#exp=lcpm
exp=expr_batch
S1PR3 <- data.frame(exp=as.numeric(exp[rownames(exp)=='S1PR3',]),
                    tissue=tissue,
                    row.names = colnames(exp))

S1PR3$tissue <- factor(S1PR3$tissue,levels = c('TALL','HSC','Thy','Tcell'))
S1PR3

library(ggplot2)
library(ggpubr)
library(ggsignif)
mycol <- c('#FF0000','#EC7D31','#00AF50','#5B9AD4')
names(mycol) <- c('TALL','HSC','Thy','Tcell')

p <- ggplot(S1PR3, aes(x=tissue, y=exp,fill=tissue)) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1)+
  geom_jitter(width = 0.1) +
  theme_classic()+labs(title="S1PR3",y='The expression level',x='tissue')+
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.title  = element_text(size=18,face="plain",color="black"),
        legend.position="none")+
  geom_signif(comparisons = list(c("TALL", "HSC"),
                                 c("TALL", "Thy"),
                                 c("TALL", "Tcell")),
              y_position=c(5.2, 5.5, 5.8,6.1),
              map_signif_level=TRUE)+
  scale_fill_manual(values =mycol)
p
table(tissue)
ggsave('SY/res/S1PR3.pdf',p,width = 8,height = 8)

