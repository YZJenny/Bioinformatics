### IC50 的散点图
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

cellline=c('PEER','CCRF_CEM','P12',"MOLT3","MOLT4",'CUTTL1','JURKAT')

df <- data.frame(cellline=cellline,
                 IC50=log(c(80.92,66.73,109.1,392.7,903.1,108.5,109.4),10))
df$group=c('red','red','black','black','black','black','black')

df$cellline=factor(df$cellline,levels = cellline)


p1 <- ggplot(df,aes(x=cellline,y=IC50)) + 
  geom_point()+theme_bw()+
  labs(x=NULL,y=NULL,fill=NULL)+  
  coord_cartesian(ylim = c(60,150))+
  theme(axis.text = element_text(color='black'),
        panel.grid.major =element_blank (), 
        panel.grid.minor = element_blank (), 
        panel.background = element_blank ())
#画上面
p2 <- ggplot(df,aes(x=cellline,y=IC50,fill=cellline)) +
  geom_point()+theme_bw()+
  labs(x=NULL,y=NULL,fill=NULL) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text = element_text(color='black')) +   
  coord_cartesian(ylim = c(350,950)) +
  scale_y_continuous(breaks = c(350,950,50))+
  theme(panel.grid.major =element_blank (), 
        panel.grid.minor = element_blank (), 
        panel.background = element_blank ())
#拼起来
p.merge <- ggarrange(p2,p1,heights=c(1/3, 2/3),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
ggsave('/mdshare/node9/yanzijun/CRU/TALL_FM/FIG_revision/IC50_TUBA1A.pdf',width = 5,height = 1.5)


