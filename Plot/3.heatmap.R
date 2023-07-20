
##################
## heatmap:subtype
##################
TAL1 <- new_info$sampleID[new_info$subtype=='TAL1']
TAL2 <- new_info$sampleID[new_info$subtype=='TAL2']
HOXA <- new_info$sampleID[new_info$subtype=='HOXA']
LMO1_2 <- new_info$sampleID[new_info$subtype=='LMO1/2']
LMO2_LYL1 <- new_info$sampleID[new_info$subtype=='LMO2_LYL1']
NKX2_1 <- new_info$sampleID[new_info$subtype=='NKX2_1']
TLX1 <- new_info$sampleID[new_info$subtype=='TLX1']
TLX3 <- new_info$sampleID[new_info$subtype=='TLX3']
Unknown <- new_info$sampleID[new_info$subtype=='Unknown']

sample_order <- c(TAL1,TAL2,HOXA,LMO1_2,LMO2_LYL1,NKX2_1,TLX1,TLX3,Unknown )
plot.exp <- dplyr::select(as.data.frame(t(candi.EXP)),sample_order) #列匹配
print(all(colnames(plot.exp)==sample_order))


library(RColorBrewer)
library(ggplot2)
col2 =brewer.pal(9,"Set1")[1:9]

group_df = data.frame(
  Subtype=rep(c('TAL1','TAL2','HOXA','LMO1_2','LMO2_LYL1','NKX2_1','TLX1','TLX3','Unknown' ), 
              c(length(TAL1),length(TAL2),length(HOXA),
                length(LMO1_2),length(LMO2_LYL1),length(NKX2_1),length(TLX1),length(TLX3),length(Unknown))))
rownames(group_df) <- colnames(plot.exp)
ann_colors = list(
  Subtype=c(TAL1=col2[1],TAL2=col2[2],HOXA=col2[3],LMO1_2=col2[4],LMO2_LYL1=col2[5],
            NKX2_1=col2[6],TLX1=col2[7],TLX3=col2[8],Unknown=col2[9]))


p.TARGET <- pheatmap::pheatmap(plot.exp,scale = 'row',cluster_cols = F,cluster_rows = T, 
                               annotation_col = group_df,show_colnames = F,
                               treeheight_row=0,treeheight_col=0,
                               annotation_colors = ann_colors,border_color = "white",
                               color =colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-4,4,by = 0.1))),
                               gaps_col = cumsum(
                                 c(length(TAL1),length(TAL2),length(HOXA),
                                   length(LMO1_2),length(LMO2_LYL1),length(NKX2_1),length(TLX1),length(TLX3),length(Unknown))
                               ),
                               breaks = seq(-4,4,by = 0.1),legend_breaks = seq(-4,4,2))
ggsave(paste('/mdshare/node9/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/TARGETsubtype_',cutoff,'_heatmap.pdf',sep=''),p.TARGET,width = 10,height = 4)
