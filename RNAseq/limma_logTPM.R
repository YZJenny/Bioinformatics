library(limma)
library(edgeR)

Tumor <- readRDS('/Users/yzj/Desktop/TCGA_assemble/BRCA/BRCA-TPM.mat.case.gn.RDS')
Normal <- readRDS('/Users/yzj/Desktop/TCGA_assemble/BRCA/BRCA-TPM.mat.ctl.gn.RDS')

pair_sample <- intersect(colnames(Tumor),colnames(Normal))
Tumor <- Tumor[,which(colnames(Tumor) %in% pair_sample)]
Normal <- Normal[,which(colnames(Normal) %in% pair_sample)]
colnames(Tumor) <- paste('T_',colnames(Tumor),sep='')
colnames(Normal) <- paste('N_',colnames(Normal),sep='')
EXP <- cbind(Tumor,Normal)

group_list <- as.factor(c(rep('T',ncol(Tumor)),rep('N',ncol(Normal))))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(T-N, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
res<-output[,c(1,4,5)]
up<-res[which(res[,3]<0.05 & res[,1] > 1),]
dim(up)
saveRDS(up,'/Users/yzj/Desktop/TCGA_assemble/BRCA/UP_gene.RDS')

