library(ggplot2)
library(pheatmap)
library(dplyr)

exp<- read.table('Desktop/0515/DATA/combined_fpkm.txt',header=T, sep='\t',stringsAsFactor=F)[,-c(1,7,8,9,10)]
exp$sum <- rowSums(exp[,2:9])#delete line which is all 0
exp <- exp[which(exp$sum != 0 & exp$sum > 1e-100 ),][,-10]
res <- read.table('Desktop/0515/DATA/gene_exp.txt_new',header=T, sep='\t',stringsAsFactor=F)[,c(3,13,15)]

height <- function(df){
  if(nrow(df) >= 1 & nrow(df) <= 20 ){
    height=7
  } else if(nrow(df) > 20 & nrow(df) <= 100 ) {
    height=20
  } else if(nrow(df) > 100 & nrow(df) <= 500 ) {
    height=90
  } else if (nrow(df) > 500 ) {
    height=140
  }
  return(height)
} 

#2.DEG
fdr=c(0.05,0.1)
fc_cutoff=c(1,1.3,1.5,2)
for(i in 1:length(fdr)){
  for(j in 1:length(fc_cutoff)){
    new_res <- mutate(res, sig=ifelse((res$q_value < fdr[i] & res$log2FC_new > log(fc_cutoff[j],2))| (res$q_value < fdr[i] & res$log2FC_new < -log(fc_cutoff[j],2)) ,ifelse(res$log2FC_new > log(fc_cutoff[j],2),'UP','DOWN'),'no'))
    gene_fc <- new_res[which(new_res$sig!='no'),c(1,3)]
    dim(gene_fc)
    sig_gene=gene_fc$gene
    up_gene=new_res[which(new_res$sig=='UP'),1]
    dn_gene=new_res[which(new_res$sig=='DOWN'),1]
    write.table(sig_gene,file=paste('Desktop/0515/DATA/DEG/fdr_',fdr[i],'/all_',fc_cutoff[j],'.txt',sep=''),sep="\t", quote=F, row.names=F, col.names=F)
    write.table(up_gene,paste('Desktop/0515/DATA/DEG/fdr_',fdr[i],'/up_',fc_cutoff[j],'.txt',sep=''),sep="\t", quote=F, row.names=F, col.names=F)
    write.table(dn_gene,paste('Desktop/0515/DATA/DEG/fdr_',fdr[i],'/dn_',fc_cutoff[j],'.txt',sep=''),sep="\t", quote=F, row.names=F, col.names=F)
    
    df <- exp[!duplicated(exp[,1]),]
    df <- merge(df,gene_fc,by.x = 'Symbol',by.y = 'gene',all=F)
    df <- df[!duplicated(df[,1]),]
    df <- df[order(df$log2FC_new,decreasing = T),]
    rownames(df) <- df$Symbol
    df <- df[,-c(1,10)]
    df <- as.matrix(df)
    h=height(df)
    pdf(paste('Desktop/0515/0.Data/2.HEATMAP/2.DEG/fdr_',fdr[i],'/fc_',fc_cutoff[j],'.pdf',sep=''),6,h)
    tmp=df
    #tmp=apply((tmp+1),2,log,2)
    stmp=t(apply(tmp,1,scale))
    rownames(stmp)=rownames(df)
    colnames(stmp)=colnames(df)
    this_diff=apply(stmp[,c(1:4)],1,sum)-apply(stmp[,c(5:8)],1,sum)
    stmp=stmp[order(this_diff),]
    pheatmap(stmp,cexRow=0.5,fontsize=6,fontsize_row = 4,fontsize_col=6,scale = "none",cluster_col = F, cluster_row = F,treeheight_row=0, treeheight_col=0,border_color=NA,color = colorRampPalette(colors = c("yellow2","blue2"))(100))
     dev.off()
  }
}
