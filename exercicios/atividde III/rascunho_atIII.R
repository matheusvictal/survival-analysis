# Limpeza do ambiente de trabalho
rm(list=ls(all=TRUE))

#Bibliotecas 
library(survminer)
library(survival)
library(ggfortify)
library(msm)


###########################################################################

# Ex 1-)

# a)


rm(list=ls(all=TRUE))

library(survminer)
library(survival)
library(ggfortify)
library(msm)


t <- c(8,11,19,24,28,33,36,38,44,96,124,130,250,250,250,7,8,10,12,13,14,19,23,25,
       26,27,31,31,49,59,64,87,89,107,117,119,
       130,148,153,156,159,191,222,rep(250,16))
d <- c(1,1,1,0,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,1,1,0,1,0,0,1,1,1,1,1,1,
       1,1,1,1,1,1,rep(0,16))

grupo <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
           2,2,2,2,2,2,2,2,2)


fit <- survfit(Surv(t,d) ~ grupo)

autoplot(fit, conf.int = F)


t1 <- c(8,11,19,24,28,33,36,38,44,96,124,130,250,250,250)
t2 <- c(7,8,10,12,13,14,19,23,25,26,27,31,31,49,59,64,87,89,107,117,119,
       130,148,153,156,159,191,222,rep(250,16))

d1 <- c(1,1,1,0,1,1,0,1,1,1,1,1,1,0,0)
d2 <- c(1,1,1,1,1,0,1,1,0,1,1,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,rep(0,16))


fit1 <- survfit(Surv(t1,d1) ~ 1)
fit2 <- survfit(Surv(t2,d2) ~ 1)


plot(log(fit1$time),log(-log(fit1$surv)) ,xlab="log (t)",
     ylab="log(-log(S(t)))", 
     ylim=c(-4,2),lty=1, lwd=2.5,pch=1,col=1)

abline(lm(log(-log(fit1$surv))~log(fit1$time)), col = 1)

points(log(fit2$time),log(-log(fit2$surv)) ,xlab="log (t)",
       ylab="log(-log(S(t)))", 
       lty=2, lwd=2.5,pch=2,col=2)

abline(lm(log(-log(fit2$surv))~log(fit2$time)), col = 2)

legend("topleft",c("Compath","Zena"), bty="n", pch=1:3,col=1:3,lty=1:3,lwd=2.5)




# b)

fitCox <- coxph(Surv(t,d) ~ grupo)
fitCox

# c)

residuos <- resid(fitCox, type = "martingale")

m <- residuos

e <- d - m # residuos de Cox-Snell

cox_snell_sum <- summary(survfit(Surv(e,d)~1))


plot(y = cox_snell_sum$surv, x = cox_snell_sum$time, ylab="S(e)",
     xlab="Residuo-Cox-Snell",lty=1,lwd=2,main=" ", type = 's')

curve(exp(-x),0,max(cox_snell_sum$time), add=T,lwd=2, lty=2, col = "blue")


# d) 





# e)

fitWei <- survreg(Surv(t,d)~grupo, dist = "weibull")
fitWei




# f)

fitCox$loglik[2]
fitWei$loglik[2]

extractAIC(fitCox) # A
extractAIC(fitWei)

BIC(fitCox) # Bayesian Information Criterion
BIC(fitWei)



# Ex 2-)

#### Testes Acelerados ####


# Para 75 kV:

## Tempos
t1 <- c(68.85, 70, 76.75, 108, 110.29, rep(120,5))
## Censuras
c1 <- c(rep(1,5),rep(0,5))
# Fator de estresse
x1 <- c(rep(28,10))


# Para 30 kV

## Tempos
t2 <- c(32.76, 35.66, 35.76, 39.85, 40, 40.25, 47.05, 54, 72, 81)
## Censuras
c2 <- rep(1,10)
## Fator de estresse
x2 <- c(rep(30,10))


# Para 32 kV

## Tempos
t3 <- c(0.4, 0.69, 0.7, 2.75, 3.75, 3.91, 4.25, 5.75, 12, 15.93)
## Censuras
c3 <- rep(1,10)
## Fator de estresse
x3 <- c(rep(32,10))


# Observação do comportamento das curvas de sobrevivência por ajustes KM

KM_1=survfit(Surv(t1, c1) ~ 1, se.fit = FALSE)
KM_2=survfit(Surv(t2, c2) ~ 1, se.fit = FALSE)
KM_3=survfit(Surv(t3, c3) ~ 1, se.fit = FALSE)





# Ex 3-)






