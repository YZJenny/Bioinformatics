#GEO 上使用 Agilent 双色 芯片的数据集的表达谱，常常给出的是两种基因的表达强度比值，并取 log2
#双色芯片处理产生的文件，是用Genepix软件得到的
rm(list=ls())
BiocManager::install('GEOquery')
library(GEOquery)
library(limma)

############
##### 1.获取基因表达谱和meta信息
############
## 1.1 获取meta信息
gset <- getGEO(GEO = 'GSE37924',getGPL = F)
meta <- pData(phenoData(gset[[1]]))
dim(meta)
patientID=apply(as.matrix(meta$title),1,function(x) unlist(strsplit(x,' '))[1])
table(patientID)

## 1.2 构造targets文件，第一列gpr的文件名，另两列是Cy3和Cy5,表示在每个gpr文件中Cy3和Cy5两种染料标记了哪组样本
setwd('/mdshare/node9/yanzijun/Extra/230223_GSE37924/GSE37924_RAW/')
dir()

FileNames <- paste(meta$geo_accession,
              paste('L',apply(as.matrix(meta$title),1,function(x) unlist(strsplit(x,'_'))[2]),sep=''),
              'res2.gpr.gz',
              sep='_')
Cy3 <- meta$`phenotype:ch1`
Cy5 <- meta$`phenotype:ch2`

targets <- data.frame(FileNames=FileNames,Cy3=Cy3,Cy5=Cy5)
f <- function(x) as.numeric(x$Flags > -75)
RG <- read.maimages(targets, source="genepix", wt.fun=f)

## 1.3 背景校正和标准化
RGne <- backgroundCorrect(RG, method="normexp", offset=25)
MA <- normalizeWithinArrays(RGne)
Expr <- as.data.frame(MA$M,row.names = MA$genes$Name)
colnames(Expr) <- paste(colnames(Expr),'.gz',sep='')

############
##### 2. 求差异基因
############
## 2.1 区分vasospasm和no vasospasm
targets$pheno <- paste(meta$`phenotype:ch1`,meta$`phenotype:ch2`,sep=':')
type <- c()
for(i in 1:length(targets$pheno)){
  item=targets$pheno[i]
  if(item=='various:no vasospasm'|item=="no vasospasm:various"){
    type[i] <- 'noVasospasm'
  }else if(item=='various:vasospasm'|item=="vasospasm:various"){
    type[i] <- 'Vasospasm'
  }
}
targets$type=type

## 2.2 用samr做差异分析(Significance Analysis of Microarrays)
data=list(x=as.matrix(Expr),y=as.numeric(as.factor(targets$type)), 
          geneid=as.character(1:nrow(Expr)),
          genenames=rownames(Expr), 
          logged2=TRUE
)
samr.obj<-samr(data, resp.type="Two class unpaired", nperms=500)
pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
FC=samr.obj$foldchange
samr.res <- data.frame(pvalue=pv,FC=FC)
samr.res <- samr.res[order(samr.res$pvalue,decreasing = F),]

genes <- c('NAMPT','NRG1','HMGCL','HTRA1','PPP2R5C')
print(samr.res[rownames(samr.res) %in% genes,])


###
group_list <- factor(type,levels = c('noVasospasm','Vasospasm') )
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(Expr)

contrast.matrix <- makeContrasts(Vasospasm-noVasospasm, levels=design)
fit <- lmFit(Expr, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
upDEG<-output[which(output$P.Value<0.05 & output$logFC > 0),]
dim(upDEG)
dnDEG<-output[which(output$P.Value<0.05 & output$logFC < 0),]
dim(dnDEG)

## 2.3 保存结果
write.csv(Expr,'../Expr.csv',row.names = T)
write.csv(targets,'../Metainfo.csv',row.names = F)
write.csv(samr.res,'../Diff_res.csv',row.names = T)


### 验证FC方向
va <- Expr[,colnames(Expr) %in% targets$FileNames[targets$type=='Vasospasm']]
noVa <- Expr[,colnames(Expr) %in% targets$FileNames[targets$type=='noVasospasm']]

genes <- rownames(Expr)[1:4]
for(i in 1:4){
  print(genes[i])
  #相减
  FC <- mean(as.numeric(va[i,]))-mean(as.numeric(noVa[i,]))
  print(FC)
}
output[rownames(output) %in% genes,]
