library(dplyr)

counts1 <- read.table('/local/yanzijun/SY/data/GSE141140_all.counts.exp.txt',header = TRUE)
gene.len<- read.table("CRU/TALL/data/GSE33470/All_hg19gene_len.txt",header = TRUE) #基因长度

counts1<-left_join(counts1,gene.len,by="Gene")#根据基因那列进行合并
counts1 <- na.omit(counts1)#删除错误值行

rownames(counts1)<-counts1[,1]
counts1<-counts1[,-1]
head(counts1)#最后一列Length是基因长度

#TPM计算
kb <- counts1$Length / 1000
countdata <- counts1[,1:(ncol(counts1)-1)]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
tpm <- as.matrix(tpm)
log.tpm1 <- log(tpm+1,2)
write.table(tpm,file="/local/yanzijun/SY/data/GSE141140_all.ltpm.exp.txt",sep="\t",quote=F)
