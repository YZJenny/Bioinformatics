#java -cp gsea-3.0.jar -Xmx4096m xtools.gsea.Gsea -out WT_KO/ -rnd_seed 13579 -set_min 5 -set_max 1000 -res WT_KO/combined_fpkm.txt_WT_KO.human.pure.des.txt -cls WT_KO/WT_KO.cls -gmx WT_KO/pathway.gmt -plot_top_x 2000 -permute gene_set -collapse false
library(ggplot2)
library(forcats)

A<-read.table(file="Desktop/0507/Apoptosis/DATA/2_GSEA/gsea_res/WT.xls",sep="\t",head=T,row.names=1,fill=T)
B<-read.table(file="Desktop/0507/Apoptosis/DATA/2_GSEA/gsea_res/KO.xls",sep="\t",head=T,row.names=1,fill=T)

gsea_dot <- rbind(A,B)
dim(gsea_dot)
summary(gsea_dot$NES)
p<- ggplot(gsea_dot, aes(x = NES, y = fct_reorder(GS.br..follow.link.to.MSigDB, NES))) + 
  geom_point(aes(color = FDR.q.val,size = SIZE)) +
  coord_cartesian(xlim=c(-2,2))+ 
  scale_x_continuous(breaks=c(-2,-1,0,1,2))+
  scale_colour_gradientn(limits=c(0, 1), colours=rainbow(6)) +
  theme_bw(base_size = 14) +
  ylab(NULL) +scale_fill_brewer(palette = 'Accent')+theme(panel.background=element_rect(fill='grey95'),panel.border = element_blank())+
  geom_vline(xintercept = 0,size=1)+
  theme(axis.text=element_text(size=3),axis.text.x = element_text(size = 9))+
  geom_segment(aes(x=-2,y=0,xend=2,yend=0))

pdf(file="Desktop/0507/Apoptosis/DATA/2_GSEA/Figure/gsea_dotplot.pdf",8,10) 
p
dev.off() 
