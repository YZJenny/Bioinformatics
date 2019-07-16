
####单因素Cox
sur_df <- as.data.frame(cbind(surtime,surevent,t(exp)))

covariates <- rownames(exp)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(surtime, surevent)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = sur_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=4)
                         wald.test<-signif(x$wald["test"], digits=4)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)



####多因素Cox
sig_pvalue<- pvalue[which(pvalue<0.05)]
sur_df <- as.data.frame(cbind(surtime,surevent,t(exp[which(rownames(exp) %in% names(sig_pvalue)),])))
f <- as.formula(paste("Surv(surtime, surevent) ~ ", paste(names(sig_pvalue), collapse = "+")))
res.cox <- coxph(f, data =  sur_df)
res.cox
