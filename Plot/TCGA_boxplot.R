## LUAD stage ia/ib/iia/iib差异
data <- read.table('/mdshare/node9/yzj/publicData/TCGAbiolinks/LUAD/LUAD-FPKM.mat.txt.case',
                   header = T,row.names = 1)
meta <- read.table('/mdshare/node9/yzj/publicData/TCGAbiolinks/LUAD/clinical.tsv',header = T,
                   fill = TRUE,sep='\t')

tmp <- apply(as.matrix(colnames(data)),1,function(x) 
  paste(unlist(strsplit(x,'\\.'))[1:2],collapse ='-')
  )
newcol <- apply(as.matrix(tmp),1,function(x) gsub('X','TCGA-',x))
colnames(data) <- newcol

sub.meta <- meta[meta$tumor_stage %in% c('stage ia','stage ib','stage iia','stage iib'),]

olsample <- intersect(colnames(data),sub.meta$submitter_id)

ol.meta <- sub.meta[sub.meta$submitter_id %in% olsample,]
ol.data <- select(data,ol.meta$submitter_id)
print(all(ol.meta$submitter_id == colnames(ol.data)))

## ENSG2genename
ENSG2name <- read.table('/mdshare/node9/yzj/publicData/TCGAbiolinks/gencode2ensg2gene.txt')
colnames(ENSG2name) <- c('ENSG','Gene')
length(intersect(ENSG2name$ENSG,rownames(data)))
LPP2.ENSG <- ENSG2name$ENSG[ENSG2name$Gene %in% c('LPP2','PPAP2C','PLPP2')]

##LPP2/PPAP2C/PLPP2
df <- data.frame(sampleID=ol.meta$submitter_id,
                 LPP2=as.numeric(ol.data[rownames(ol.data)==LPP2.ENSG,]),
                 Stage=ol.meta$tumor_stage)
df$LPP2 <- log(df$LPP2+1,2)

library(ggplot2)
library(ggpubr)
library(ggsignif)
mycol <- c('#FF0000','#EC7D31','#00AF50','#5B9AD4')
names(mycol) <- c('stage ia','stage ib','stage iia','stage iib')

p <- ggplot(df, aes(x=Stage, y=LPP2,fill=Stage)) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1)+
  geom_jitter(width = 0.0001) +
  theme_classic()+labs(title="LPP2",y='The expression level',x='Tumor stage')+
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        axis.text = element_text(size=18,face="plain",color="black"),
        axis.title  = element_text(size=18,face="plain",color="black"),
        legend.position="none")+
  geom_signif(comparisons = list(c("stage ia", "stage ib"),
                                 c("stage ia", "stage iia"),
                                 c("stage ia", "stage iib")),
              y_position=c(5.2, 5.5, 5.8,6.1),
              map_signif_level=TRUE)+
  scale_fill_manual(values =mycol)
p

library(ggsci)
p1<- ggplot(df,aes(Stage,LPP2,fill=Stage)) + 
  geom_boxplot()+
  scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  #stat_compare_means(method = "anova")+
  stat_compare_means()+theme_bw()

p2 <- ggplot(df,aes(Stage,LPP2,fill=Stage)) + 
  geom_boxplot()+
  scale_fill_jco()+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  geom_signif(comparisons = list(c("stage ia", "stage ib"),c("stage ia", "stage iia"),
                                 c("stage ia", "stage iib")),
              y_position=c(5.5, 5.8, 6.1,6.4))+theme_bw()
pdf('/mdshare/node9/yanzijun/Extra/230224_LUAD/LPP2_Stage.pdf',width = 5,height = 4)
print(p1)
print(p2)
dev.off()

