library(ggplot2)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)

exp_NG2017 <- read.table('CRU/TALL/data/NG_2017/inputExp.txt',header = T,sep='\t')
EFS <- read.table('CRU/TALL/data/NG_2017/EFS_TAL.txt',header = T,
                  row.names = 1,sep='\t')
EFS_TAL <- EFS[EFS$surtype=='TAL',]
dim(EFS_TAL)
expEFS_TAL <- dplyr::select(exp_NG2017,rownames(EFS_TAL))
dim(expEFS_TAL)
print(all(colnames(expEFS_TAL)==rownames(EFS_TAL)))


OS <- read.table('CRU/TALL/data/NG_2017/OS_TAL.txt',header = T,
                 row.names = 1,sep='\t')
OS_TAL <- OS[OS$surtype=='TAL',]
dim(OS_TAL)
expOS_TAL <- dplyr::select(exp_NG2017,rownames(OS_TAL))
dim(expOS_TAL)
print(all(colnames(expOS_TAL)==rownames(OS_TAL)))


EFS_plots <- list()
OS_plots <- list()

for(i in 1:length(Genes)){
  #### EFS
  exp_TAL=expEFS_TAL
  sur_df <- EFS_TAL
  
  RS=as.numeric(exp_TAL[rownames(exp_TAL)==Genes[i],])
  cut_off=mean(RS)
  surtype=rep('high exp',length(RS))
  surtype[which(RS < cut_off)]='low exp'
  
  sur_df$surtype <- surtype
  fit <- survfit(Surv(time =time, event =status) ~ surtype,data=sur_df)
  print(surv_pvalue(fit)$pval)
  p <- ggsurvplot(fit, pval = TRUE,ggtheme = theme_survminer(base_size = 10),title=Genes[i])
  EFS_plots[[i]] <- p
  
  #### OS
  exp_TAL=expOS_TAL
  sur_df <- OS_TAL
  
  RS=as.numeric(exp_TAL[rownames(exp_TAL)==Genes[i],])
  cut_off=mean(RS)
  surtype=rep('high exp',length(RS))
  surtype[which(RS < cut_off)]='low exp'
  
  sur_df$surtype <- surtype
  fit <- survfit(Surv(time =time, event =status) ~ surtype,data=sur_df)
  print(surv_pvalue(fit)$pval)
  p <- ggsurvplot(fit, pval = TRUE,ggtheme = theme_survminer(base_size = 10),title=Genes[i])
  OS_plots[[i]] <- p
}
# 将多个图合并一起
p <- arrange_ggsurvplots(EFS_plots, print = TRUE) #定义行数和列数
ggsave(filename=paste('/local/yanzijun/SY/res/EFS_bloodGenes.pdf'),p)

p <- arrange_ggsurvplots(OS_plots, print = TRUE) #定义行数和列数
ggsave(filename=paste('/local/yanzijun/SY/res/OS_bloodGenes.pdf'),p)
