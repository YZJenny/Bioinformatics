## 分组棒棒图，参考https://zhuanlan.zhihu.com/p/480622216
library(forcats)
library(ggstance)
library(ggplot2)
library(RColorBrewer)
my.col <- brewer.pal(12,"Set3")[-2]

GSEA.df <- readRDS("Hallmarks/gsea_allcancer.rds")

targets <- c("HALLMARK_ANGIOGENESIS", "HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_G2M_CHECKPOINT",
             "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_MTORC1_SIGNALING")

for(i in 1:length(targets)){
  t <- targets[i]
  print(t)
  title <- str_to_title(gsub("HALLMARK_", "", t))
  
  data <- GSEA.df[GSEA.df$ID==t,]
  data$descrip <- paste(data$pathway,data$cancer,sep='_')
  list <- sort(data$descrip,decreasing = TRUE)
  list <- factor(list,levels = list)
  data$logPvalue <- -log(data$p.adjust, 10)
  
  
  p <- ggplot(data,aes(x=NES,y=factor(descrip,levels = list)))+
    geom_point(aes(size=logPvalue,color=pathway))+
    geom_segment(aes(x=0,xend=NES,
                     y=descrip,yend=descrip,
                     color=pathway),
                 cex=1)+
    scale_y_discrete(position = 'right')+
    labs(title = title,y=NULL,color='Pathway Category',size='-log(p.adjust,10)')+
    theme_bw()+
    #scale_size(range = c(1,16),breaks = c(1.3,5,10,15))+ #自定义点大小图例的刻度范围
    scale_color_manual(values = my.col[1:7])+
    theme(plot.title = element_text(size=12,hjust=0.5), #标题居中
          plot.margin = margin(0.5,0.5,4.5,0.5,'cm'),
          legend.position = c(0.9,-0.17),
          legend.box = 'horizontal',
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))+
    guides(size=guide_legend(ncol = 4,order = 0,
                             label.position = 'bottom'),
           color=guide_legend(ncol = 2,order=1))+
    theme(axis.ticks.length.x = unit(0.05,'cm'),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=10,colour = 'black'),
          axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 10,colour = 'black'
                                     #face = 'bold'
          ))
  ggsave(paste('Hallmarks/barplot_',title,'.pdf',sep=''),p,width = 6,height = 10)
  
}
