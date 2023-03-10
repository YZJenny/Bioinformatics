## 自己ggplot画图
rm(list=ls())
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(ggplot2)

setwd('/mdshare/node9/yanzijun/Extra/230310_GOKEGG/')
geneLst <- read.table('gene.txt')[,1]
print(length(geneLst))

symbol2id=mapIds(org.Hs.eg.db,geneLst,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
print(length(id))

#GO分析#
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "all") #GO富集分析
print(dim(ego))
ego_res <- as.data.frame(ego)
write.csv(ego_res,'GO.csv',row.names = F)

#KEGG分析#
ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
ekk_res <- as.data.frame(ekk)
dim(ekk_res)
write.csv(ekk_res,'KEGG.csv')

##画图
library(RColorBrewer)
my_palette <-colorRampPalette(rev(brewer.pal(9,"YlOrRd")),alpha=TRUE)(n=399)
get_plot <- function(mat.plot){
  mat.plot$logpadjust <- -log(mat.plot$p.adjust,10)
  mat.plot$GeneRatio <- apply(as.matrix(mat.plot$GeneRatio),1,function(x) 
    round(as.numeric(unlist(strsplit(x,'/'))[1])/as.numeric(unlist(strsplit(x,'/'))[2]),2))
  mat.plot <- mat.plot[order(mat.plot$p.adjust,decreasing = F),]
  mat.plot$Description <- factor(mat.plot$Description, levels = rev(mat.plot$Description))
  
  plot <- ggplot(mat.plot,aes(x=GeneRatio,y=Description)) +
    geom_point(aes(size=Count,color=logpadjust)) +
    scale_color_gradientn('-log(p.adjust,10)', colors=my_palette) +
    theme_bw() +
    theme(axis.text=element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))+
    labs(x="Gene Ratio")
  return(plot)
}

pickterm <- read.table('画图.txt',header = F,sep='\t')
pickBP <- pickterm$V3[pickterm$V1=='BP']
pickCC <- pickterm$V3[pickterm$V1=='CC']

mat.plot <- ego_res[ego_res$Description %in% pickBP,c('Description','GeneRatio','p.adjust','Count')]
BP.plot <-get_plot(mat.plot)
BP.plot <- BP.plot+labs(title="BP Enrichment")
ggsave('BP.pdf',BP.plot,width = 9,height = 4)

mat.plot <- ego_res[ego_res$Description %in% pickCC,c('Description','GeneRatio','p.adjust','Count')]
CC.plot <-get_plot(mat.plot)
CC.plot <- CC.plot+labs(title="CC Enrichment")
ggsave('CC.pdf',CC.plot,width = 6,height = 3)

mat.plot <- ekk_res[,c('Description','GeneRatio','p.adjust','Count')]
KEGG.plot <-get_plot(mat.plot)
KEGG.plot <- KEGG.plot+labs(title="KEGG Pathway Enrichment")
ggsave('KEGG.pdf',KEGG.plot,width = 9,height = 4)
