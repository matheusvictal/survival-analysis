library("flexsurv")
library("data.table")
library("ggplot2")
library("muhaz")
library("survival")

rm(list = ls())

y=c(3,5,6,8,10,11,15,20,22,23,27,29,32,35,40,26,28,33,21,24)
delta=c(1,1,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,0)

mKM = survfit(Surv(y, delta) ~ 1)

d <- data.frame(time = mKM$time,
                n.risk = mKM$n.risk,
                n.event = mKM$n.event,
                n.censor =mKM$n.censor,
                surv = mKM$surv,
                upper = mKM$upper,
                lower = mKM$lower)
                
fit_w <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "Weibullph") 
                
fit_G <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "gamma") 
                
fit_LL <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "llogis") 
                
fit_LN <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "lognormal") 
                
fit_Go <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "gompertz") 
                
fit_GG <- flexsurvreg(Surv(y, delta) ~ 1,  dist = "gengamma") 
                

               
par(mfrow=c(2,3))
                
plot(fit_w, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")
                
title(main="Weibull")
                
plot(fit_G, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")
                
title(main="Gama")
                
plot(fit_LL, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")
                
title(main="Log-Logistic")
                
plot(fit_LN, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")
                
title(main="Log-Normal")
                  
plot(fit_Go, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")

title(main="Gompertz")

plot(fit_GG, ci=FALSE, conf.int=FALSE, ylab="Survival", xlab="Years")

title(main="Generalized Gamma")



crit=rbind(Wei=glance(fit_w),glance(fit_G),
           glance(fit_LL),glance(fit_LN),
           glance(fit_Go),glance(fit_GG))     
library(xtable)
xtable(crit)
