rm(list=ls())
library(ggplot2)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)

setwd('/mdshare/node9/yanzijun/Extra/230224_surv_cutoff/')
data <- read.csv('data.csv',row.names = 1)

##surv_cutpoint
Genes=colnames(data)[1:8]
res.cut <- surv_cutpoint(data,time = "OS.time", event = "OS",variables = Genes)
cut.df <- as.data.frame(summary(res.cut))
cut.df <- tibble::rownames_to_column(cut.df,'gene')
res.cat <- surv_categorize(res.cut)
write.csv(cut.df,'cutpoint.csv',row.names = F)
write.csv(res.cat,'category.csv',row.names = T)

clin_info=res.cat[,1:2]
group_info <- res.cat[,3:ncol(res.cat)]
pvalue <- c()
for(k in 1:length(Genes)){
  print(k)
  surtype=group_info[,k]
  clin_info$surtype <- surtype
  fit <- survfit(Surv(time =OS.time, event =OS) ~ surtype,data=clin_info)
  pvalue[k] <- surv_pvalue(fit)$pval
}
names(pvalue) <- Genes
print(length(pvalue[pvalue <0.05]))
print(pvalue)


##median/mean
tmp <- function(x){
  #survtype <- rep('Low',length(x))
  #survtype[x>quantile(x,0.5)] <- 'High'
  
  # survtype <- rep('median',length(x))
  # survtype[x>quantile(x,0.75)] <- 'High'
  # survtype[x<quantile(x,0.25)] <- 'Low'
  
  survtype <- rep('Low',length(x))
  survtype[x>mean(x)] <- 'High'
  return(survtype)
}
group_info <- as.data.frame(apply(data[,1:8],2,tmp))
clin_info=data[,9:10]

pvalue <- c()
for(k in 1:length(Genes)){
  print(k)
  surtype=group_info[,k]
  clin_info$surtype <- surtype
  fit <- survfit(Surv(time =OS.time, event =OS) ~ surtype,data=clin_info)
  pvalue[k] <- surv_pvalue(fit)$pval
}
names(pvalue) <- Genes
print(pvalue)
