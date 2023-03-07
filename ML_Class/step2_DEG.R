rm(list=ls())
setwd('/mdshare/node9/yanzijun/Extra/230228_Model/')

############
### 1. limma, EXP是log后的表达矩阵（以2为底）
############
library(limma)
data.lst <- readRDS('input_method3/data.lst')

EXP <- data.lst$train.exp
if(max(EXP) > 50){
  EXP <- log(EXP+1,2)
}

Label <- data.lst$train.label
Label$lable <- 'Metastasis'
Label$lable[Label$Group==0] <- 'noMetastasis'
Label <- Label$lable

group_list <- factor(Label,levels=c('Metastasis','noMetastasis'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(Metastasis-noMetastasis, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene') 
min(output$adj.P.Val)
write.table(output,'Fig_method3/allDiff.txt',sep='\t',row.names = F,quote = F)
sigG <- output$gene[which(output$adj.P.Val < 0.05 & abs(output$logFC) > log(1,2))]
print(length(sigG))

############  
### 2. PCA
############
rm(list=ls())
library(ggplot2)
library("FactoMineR")
library("factoextra")
setwd('/mdshare/node9/yanzijun/Extra/230228_Model/')
output <- read.table("Fig_method3/allDiff.txt",header = T,sep='\t') 
sigG <- output$gene[which(output$P.Value < 0.05 & abs(output$logFC) > log(2,2))]
print(length(sigG))

data.lst <- readRDS('input_method3/data.lst')
EXP <- data.lst$train.exp
Label <- data.lst$train.label
Label$lable <- 'Metastasis'
Label$lable[Label$Group==0] <- 'noMetastasis'
Label <- Label$lable

sigEXP <- log(EXP[sigG,]+1,2)
newcol <- apply(as.matrix(rownames(sigEXP)),1,function(x) gsub('-','',x))
rownames(sigEXP) <- newcol

df.pca <- PCA(t(sigEXP), graph = FALSE)
pca.plot <- fviz_pca_ind(df.pca,
             geom.ind = "point",
             col.ind = factor(Label) ,
             addEllipses = TRUE,
             legend.title = "Groups"
)

pdf('Fig_method3/PCA.pdf',width = 6,height = 4)
print(pca.plot)  ## DEG
dev.off()


############  
### 3. vlnplot
############
rm(list=ls())
library(ggplot2)
library(ggrepel)
library(dplyr)
data <- read.table("Fig_method3/allDiff.txt",header = T,sep='\t') 

qvalue <- 0.05
fc_cutoff <- 1.5
data <- mutate(data, 
               sig=ifelse((data$adj.P.Val < qvalue & data$logFC > log(fc_cutoff,2))| 
                            (data$adj.P.Val < qvalue & data$logFC < -log(fc_cutoff,2)) ,
                          ifelse(data$logFC > log(fc_cutoff,2),'UP','DOWN'),'no'))
data <- na.omit(data)
data$log10Qvalue=-log10(data$adj.P.Val)
table(data$sig)

UP.data <- data[data$sig=='UP',]
DN.data <- data[data$sig=='DOWN',]

UP.data.FC <- UP.data[order(UP.data$logFC,decreasing = T),]
DN.data.FC <- DN.data[order(DN.data$logFC,decreasing = F),]

tag_lst=c(head(UP.data$gene,10),head(DN.data$gene,10),
          head(UP.data.FC$gene,10),head(DN.data.FC$gene,10))
tag_gene=as.character(data$gene)
index_tag=which(!tag_gene %in% tag_lst)
tag_gene[index_tag]=""

figure <- ggplot(data,aes(x=logFC,y=log10Qvalue))+geom_point(aes(color=sig))+
  scale_color_manual(values = c("#999999","#FF0000"))+
  #scale_color_manual(values = c("#999999"))+
  labs(x="Log2(Fold Change)",y = "-Log10(adj.P value)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=10),panel.grid.major=element_line(colour=NA))+
  theme_bw()+
  geom_vline(xintercept=log(fc_cutoff,2),lty=3,col="black",lwd=0.5)+
  geom_vline(xintercept=-log(fc_cutoff,2),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(qvalue),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
figure
ggsave(plot=figure,filename='Fig_method3/volcano.pdf',width = 8,height =8 )

