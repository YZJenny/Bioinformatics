library(ggplot2)
library(ggrepel)
library(dplyr)

data <- read.csv("allDiff.csv",header = T)
colnames(data)[1] <- 'gene'

qvalue <- 0.05
fc_cutoff <- 1
data <- mutate(data, 
               sig=ifelse((data$adj.P.Val < qvalue & data$logFC > log(fc_cutoff,2))| 
                            (data$adj.P.Val < qvalue & data$logFC < -log(fc_cutoff,2)) ,
                          ifelse(data$logFC > log(fc_cutoff,2),'UP','DOWN'),'no'))
data <- na.omit(data)
data$log10Qvalue=-log10(data$adj.P.Val)


tag_lst=paste("S1PR",1:5,sep ='')
tag_gene=as.character(data$gene)
index_tag=which(!tag_gene %in% tag_lst)
tag_gene[index_tag]=""


figure <- ggplot(data,aes(x=logFC,y=log10Qvalue))+geom_point(aes(color=sig))+
  scale_color_manual(values = c("#5757F8","#999999","#FF0000"))+
  labs(x="Log2(fold change)",y = "-Log10(Q value)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=7),panel.grid.major=element_line(colour=NA))+
  theme_bw()

p <- figure+geom_text_repel(label=tag_gene,max.overlaps = 30)+
  geom_vline(xintercept=0,lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p
ggsave(plot=p,filename = 'volcano.pdf',width = 5,height = 5)
