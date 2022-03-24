library(ggplot2)
library(pheatmap)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(rmda)
library(DCA)
get_Unires <- function(Cox_df){
  covariates <- colnames(Cox_df)[4:ncol(Cox_df)]
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Cox_df)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$sctest["pvalue"], digits=2)
                           sc.test<-signif(x$sctest["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR.confint <- paste0(HR, " (", 
                                                HR.confint.lower, "-", HR.confint.upper, ")")
                           # res<-c(beta, HR.confint, sc.test, p.value)
                           # names(res)<-c("beta", "HR (95% CI for HR)", "sc.test", "p.value")
                           res<-c(HR, HR.confint.lower, HR.confint.upper, HR.confint, p.value)
                           names(res)<-c( 'HR','lower','upper',"HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res,stringsAsFactors = F)
  res$Gene <- rownames(res)
  res$adjustP <- p.adjust(res$p.value,method = 'fdr')
  res <- as.matrix(res)
  res <- rbind(colnames(res),res)
  return(res)
}
