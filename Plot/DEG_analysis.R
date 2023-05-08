rm(list=ls())
library(limma)
library(ggplot2)
library(dplyr)
setwd('/mdshare/node9/yanzijun/Extra/230503_GSE136825/')

### step1. gtf文件提取gene_id和gene_name的对应关系
gtf <- read.table('gencode.v26.annotation.gtf',header = F,sep='\t',skip = 5)
gtf[1:3,]
gtf_gene <- gtf[gtf[,3]=="gene",]
Geneid <- NULL
name <- NULL
for (i in 1:nrow(gtf_gene)){
  temp <- unlist(strsplit(gtf_gene[i,9],";"))
  Geneid[i] <- substr(temp[1],9,nchar(temp[1]))
  name[i] <- substr(temp[3],12,nchar(temp[3]))
  print(i)
}
id2name <- as.data.frame(cbind(Geneid,name))
id2name[1:3,]

### step2. 读入表达谱，并将表达谱中gene_id替换为gene_name
data <- read.table('GSE136825_genecounts_20190903.txt',header = TRUE,sep='\t')
data[1:3,1:4]
data$Chr=data$Start=data$End=data$Strand=data$Length=NULL
data[1:3,1:4]

length(intersect(data$Geneid,id2name$Geneid))

exp <- merge(id2name,data,by='Geneid')
dim(exp)
exp <- exp[!duplicated(exp$name),]
dim(exp)
rownames(exp) <- exp$name
exp$Geneid=exp$name=NULL
exp[1:4,1:5]

### step3. 提取样本类型，构建分组矩阵
gset <- getGEO("GSE136825",destdir = './')
meta <- pData(phenoData(gset[[1]]))
meta[1:3,]

case=meta$description[grep('^PY|^TT',meta$title)]
case=paste(case,'.bam',sep='')
control=meta$description[grep('^CC|^CCY|^CCS',meta$title)]
control=paste(control,'.bam',sep='')
group <- factor(c(rep('case',length(case)),rep('control',length(control))),levels = c('case','control'))

### step4. 提取case/control样本表达谱
sub.exp <- dplyr::select(exp,c(case,control))
dim(sub.exp)

design <- model.matrix(~0+group)
rownames(design) = colnames(sub.exp)
colnames(design) <- levels(group)
head(design)

### step5. 构建进行差异分析的对象
DGElist <- DGEList(counts = sub.exp, group = group)
keep_gene <- rowSums(cpm(DGElist) > 1) >= 2
table(keep_gene)
DGElist <- DGElist[keep_gene, ,keep.lib.sizes =FALSE]

### step6. 构建线性模型
DGElist <- calcNormFactors( DGElist )
# 将count值转化成log2-counts per million (logCPM)，准备进行线性回归
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
# 对每一个基因进行线性模型构建
fit <- lmFit(v, design)
# 构建比较矩阵
cont.matrix <- makeContrasts(contrasts = c('case-control'), levels = design)
# 构建芯片数据的线性模型，计算估计的相关系数和标准差
fit2 <- contrasts.fit(fit, cont.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2 <- eBayes(fit2)

### step7. 得出差异结果
nrDEG_limma_voom = topTable(fit2, coef = 'case-control', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
#筛选出符合要求的差异基因
res<-cbind(rownames(nrDEG_limma_voom),nrDEG_limma_voom)
DEG_res <-res %>% dplyr::filter((logFC>1 | logFC < (-1)) & adj.P.Val < 0.05)
colnames(res)[1]=colnames(DEG_res)[1]="Symbol"
head(DEG_res)

write.csv(res,'limma_res.csv',row.names = F)
write.csv(DEG_res,'DEG_res.csv',row.names = F)
write.csv(sub.exp,'exp.csv',row.names = F)
saveRDS(group,'group.rds')

### step8. FigA top50 DEG 热图
library(pheatmap)
group_df = data.frame(group=group)
levels(group_df$group) <- c('control','case')
rownames(group_df) <- colnames(sub.exp)

##提取top50DEG
DEG_res.order <- DEG_res[order(DEG_res$logFC,decreasing = T),]
top25.up <- DEG_res.order$Symbol[1:25]
top25.dn <- DEG_res.order$Symbol[(nrow(DEG_res.order)-24):nrow(DEG_res.order)]
top50DEG <- c(top25.up,top25.dn)

plot.exp <- log(sub.exp+1,2)
plot.exp <- plot.exp[match(top50DEG,rownames(plot.exp)),]
dim(plot.exp)
print(all(rownames(plot.exp)==top50DEG))
## 画图
p.heatmap <- pheatmap(plot.exp,scale = 'row',cluster_cols = F,cluster_rows = F,
                      annotation_col = group_df,show_colnames = F,
                      treeheight_row=0,treeheight_col=0,
                      border_color = "white",
                      colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-4,4,by = 0.1))),
                      breaks = seq(-4,4,by = 0.1),legend_breaks = seq(-4,4,2))
p.heatmap
ggsave('FigA_heatmap.pdf',p.heatmap,width = 5,height = 7)


### step9. FigB 火山图
qvalue <- 0.05
fc_cutoff <- 1
plot.res <- mutate(res, 
               sig=ifelse((res$adj.P.Val < qvalue & res$logFC > log(fc_cutoff,2))| 
                            (res$adj.P.Val < qvalue & res$logFC < -log(fc_cutoff,2)) ,
                          ifelse(res$logFC > log(fc_cutoff,2),'UP','DOWN'),'no'))
plot.res <- na.omit(plot.res)
plot.res$log10Qvalue=-log10(plot.res$adj.P.Val)

p.volca <- ggplot(plot.res,aes(x=logFC,y=log10Qvalue))+geom_point(aes(color=sig))+
  scale_color_manual(values = c("#5757F8","#999999","#FF0000"))+
  labs(x="Log2(fold change)",y = "-Log10(Q value)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=7),panel.grid.major=element_line(colour=NA))+
  theme_bw()+
  geom_vline(xintercept=0,lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p.volca
ggsave(plot=p.volca,filename = 'FigB_volcano.pdf',width = 5,height = 5)

### step10. FigC-D top10 hubDEG的degree
## 将DEG输入string网站得出PPI network，然后导入cytoscape，利用cytohubba的MCC算法，得出top10 hub gene
top.DEG_res <- DEG_res[order(abs(DEG_res.order$logFC),decreasing = F),]
top.DEG_res <- top.DEG_res[1:500,]
write.csv(top.DEG_res,'top.DEG_res.csv',row.names = F)

## FigD degree的barplot
node <- read.csv('string_interactions_node.csv',header = T)
head(node)

top10node <- c('EXO1','CDC45','KIAA0101','ASF1B','TROAP','MELK','PTTG1','UBE2C',
               'NEK2','KIF4A')
plot.df <- node[node$name %in% top10node,c('name','Degree')]
plot.df <- plot.df[order(plot.df$Degree,decreasing = T),]
plot.df$name <- factor(plot.df$name,levels = plot.df$name)
head(plot.df)

p.degree <- ggplot(plot.df,aes(name,Degree))+
  geom_bar(stat='identity',width = 0.9,colour='black',fill='cadetblue1')+
  coord_flip() + #horizontal
  labs(y='Degree',x="")+
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y= element_blank())+ #delete backgroud
  guides(fill=FALSE)
p.degree
ggsave('FigD_Degree.pdf',p.degree,width=6,height=4)

### step11. FigE-F GO/KEGG富集分析
library(org.Hs.eg.db) 
library(clusterProfiler)

geneLst <- DEG_res$Symbol
print(length(geneLst))

symbol2id=mapIds(org.Hs.eg.db,geneLst,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
print(length(id))

#GO分析#
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP") #GO富集分析
print(dim(ego))
ego_res <- as.data.frame(ego)
write.csv(ego_res,'GO.csv',row.names = F)

#KEGG分析#
ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
ekk_res <- as.data.frame(ekk)
dim(ekk_res)
write.csv(ekk_res,'KEGG.csv')

## FigE GO画图,结果太多，选取top20的terms，自己可以修改top的数目
p.GO <- barplot(ego,showCategory = 20)
ggsave('FigE_GO.pdf',p.GO,width=6,height=8)

## FigF KEGG画图,结果太多，选取top20的terms，自己可以修改top的数目
p.KEGG <- barplot(ekk,showCategory = 20)
ggsave('FigF_KEGG.pdf',p.KEGG,width=6,height=13)
