library(ggplot2)
library(ggrepel)
library(dplyr)

data <- read.table("Desktop/0515/DATA/gene_exp.txt",header = T)
data$log2FC_new <- log((data$value_2+0.01)/(data$value_1+0.01),2)#log2FC_new:value+1
data <- mutate(data, significant_new=ifelse(data$q_value < 0.1 & data$log2.fold_change.!=0,'yes','no'))#significant_new:q<0.1
print(head(data))
write.table(data,"Desktop/0515/DATA/gene_exp.txt_new",sep="\t", quote=F, row.names=F, col.names=T)

qvalue <- 0.1 #0.1/0.05
fc_cutoff <- c(1,1.3,1.5,2)
for(i in 1:length(fc_cutoff)){
  data <- mutate(data, sig=ifelse((data$q_value < qvalue & data$log2.fold_change. > log(fc_cutoff[i],2))| (data$q_value < qvalue & data$log2.fold_change. < -log(fc_cutoff[i],2)) ,ifelse(data$log2.fold_change. > log(fc_cutoff[i],2),'UP','DOWN'),'no'))
  sig=data$sig
  filter_gene=as.character(data$gene)
  filter_gene_log2FC=data$log2FC_new
  filter_gene_logP=-log10(data$p_value)
  filter_gene_qvalue=data$q_value
  filter_data=data.frame(filter_gene,filter_gene_log2FC,filter_gene_logP,filter_gene_qvalue)
  boolean_tag=(filter_gene_log2FC>log(fc_cutoff[i],2) | filter_gene_log2FC< -log(fc_cutoff[i],2)) & filter_gene_qvalue < qvalue
  if(qvalue==0.05){
    tag_lst=c("Myc","Dusp1","Fos","Jun","Pim1","Mtor")
  }else if(qvalue==0.1){
    tag_lst=c("Myc","Ccnd1","Dusp1","Fos","Jun","Pim1","Mtor")
  }
  index_tag=which(!filter_gene %in% tag_lst)
  filter_gene[index_tag]=""
  #filter_gene[which(boolean_tag==F)]=""#only tag sig gene

  figure <- ggplot(filter_data,aes(x=filter_gene_log2FC,y=filter_gene_logP))+geom_point(aes(color=sig))+
    scale_color_manual(values = c("#377EB8","#999999","#E41A1C"))+
    labs(x="Log2(fold change)",subtitle = "-Log10(p value)")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text=element_text(size=7),panel.grid.major=element_line(colour=NA))+theme_bw()+
    theme(panel.grid=element_blank(),panel.border=element_blank(),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    theme(plot.subtitle = element_text(hjust = 0.5))+
    geom_segment(aes(x=-15,y=-0.02,xend=15,yend=-0.02))+
    geom_segment(aes(x=-15,y=-0.02,xend=-15,yend=-0.1))+
    geom_segment(aes(x=-10,y=-0.02,xend=-10,yend=-0.1))+
    geom_segment(aes(x=-5,y=-0.02,xend=-5,yend=-0.1))+
    geom_segment(aes(x=0,y=-0.02,xend=0,yend=-0.1))+
    geom_segment(aes(x=5,y=-0.02,xend=5,yend=-0.1))+
    geom_segment(aes(x=10,y=-0.02,xend=10,yend=-0.1))+
    geom_segment(aes(x=15,y=-0.02,xend=15,yend=-0.1))+
    geom_segment(aes(x=0,y=-0.02,xend=0,yend=5))+
    geom_segment(aes(x=-0.3,y=1,xend=0,yend=1))+
    geom_segment(aes(x=-0.3,y=2,xend=0,yend=2))+
    geom_segment(aes(x=-0.3,y=3,xend=0,yend=3))+
    geom_segment(aes(x=-0.3,y=4,xend=0,yend=4))+
    geom_segment(aes(x=-0.3,y=5,xend=0,yend=5))+
    annotate(geom = 'text',x=-0.7,y=5,label='5')+
    annotate(geom = 'text',x=-0.7,y=4,label='4')+
    annotate(geom = 'text',x=-0.7,y=3,label='3')+
    annotate(geom = 'text',x=-0.7,y=2,label='2')+
    annotate(geom = 'text',x=-0.7,y=1,label='1')+
    annotate(geom = 'text',x=-15,y=-0.2,label='-15')+
    annotate(geom = 'text',x=-10,y=-0.2,label='-10')+
    annotate(geom = 'text',x=-5,y=-0.2,label='-5')+
    annotate(geom = 'text',x=0,y=-0.2,label='0')+
    annotate(geom = 'text',x=5,y=-0.2,label='5')+
    annotate(geom = 'text',x=10,y=-0.2,label='10')+
    annotate(geom = 'text',x=15,y=-0.2,label='15')
  p <- figure+geom_text_repel(label=filter_gene)
  ggsave(plot=p,filename = paste('Desktop/0515/0.Data/1.VOLCANO/qvalue_',qvalue,'/',fc_cutoff[i],'_withGene','.pdf',sep=''))
  ggsave(plot = figure,filename = paste('Desktop/0515/0.Data/1.VOLCANO/qvalue_',qvalue,'/',fc_cutoff[i],'_withoutGene','.pdf',sep=''))
   
}
