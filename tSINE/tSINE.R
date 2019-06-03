TAG=as.character(read.table('Desktop/GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.TAG.txt',header=F)[,1])
#a=read.table('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.txt',row.names=1,header=T,sep='\t',check.names = F)
a=readRDS('Desktop/GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.RDS')

library(Seurat)
library(dplyr)

pbmc <- CreateSeuratObject(raw.data = a, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc@meta.data$tag=TAG
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot=F,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc)


PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)

PCElbowPlot(object = pbmc)

PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, perplexity=5, dims.use = PCUSE, do.fast = TRUE)

TSNEPlot(object = pbmc,group.by='tag',pt.size=4)
PCAPlot(object = pbmc,group.by='tag',pt.size=4)

################################

SURTAG=read.table('Desktop/SURTAG.txt',sep='\t',header=F)

PC1=pbmc@dr$pca@cell.embeddings[,1]
PC1_GENE=pbmc@dr$pca@gene.loadings[,1]
GBM=which(TAG=='GBMx')
LGG=which(TAG=='LGGx')


D=which(SURTAG[,1]==1)
TAG_D=TAG[D]
COL=rep('red',length(D))
COL[which(TAG_D=='GBMx')]='indianred3'
COL[which(TAG_D=='LGGx')]='green3'

used=which(!is.na(SURTAG[,1]))
TAG_U=TAG[used]
COL_U=rep('red',length(used))
COL_U[which(TAG_U=='GBMx')]='indianred3'
COL_U[which(TAG_U=='LGGx')]='green3'

par(mfrow=c(1,2))
plot(PC1[D],SURTAG[D,2],pch=16,col=COL, xlab='PC1',ylab='OS (days)', cex=1.5)
plot(PC1[used],SURTAG[used,2],pch=16,col=COL_U, xlab='PC1',ylab='OS (days)', cex=1.5)


cor.test(PC1[D],SURTAG[D,2],method='spearman')
cor.test(PC1[used],SURTAG[used,2],method='spearman')

####################
library(survival)
library(survminer)
a=read.table('Desktop/SUR.txt',header=T,sep='\t')
used=used
score=PC1[used]
SUR=SURTAG[used,2]
SURE=SURTAG[used,1]

TYPE=rep('MED',length(used))
TYPE[which(score> quantile(score,0.5) )]='High'
TYPE[which(score<= quantile(score,0.5) )]='Low'
surtime=SUR[which(TYPE!='MED')]
surevent=SURE[(which(TYPE!='MED'))]
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)

plot(PC1[used],SURTAG[used,2],pch=16,col=COL_U, xlab='PC1',ylab='OS (days)', cex=1.5)
abline(v=quantile(score,0.5),col='red',lwd=1.5)

#############
library("cluster")
library("factoextra")
library("magrittr")
my_data <- (pbmc@dr$pca@cell.embeddings)
res.dist <- get_dist(my_data, stand = TRUE, method = "pearson")
#fviz_dist(res.dist,gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#fviz_nbclust(my_data, kmeans, method = "gap_stat")
set.seed(123)
km.res <- kmeans(my_data, 2, nstart = 25)
km.res$cluster
fviz_cluster(km.res, data = my_data,ellipse.type = "convex",palette = "jco",ggtheme = theme_minimal())

