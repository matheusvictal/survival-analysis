library(numDeriv)
library(xtable)
#####################################################
#          Weibull distribution                     #
#####################################################

dWEI4 = function (y, mu = 1, sigma = 0, log = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(y < 0))
    stop(paste("y must be greater than 0 ", "\n", ""))
  sigma2 = exp(-sigma / mu)
  fy = dweibull(y, shape = mu, scale = sigma2, log = log)
  fy
}

pWEI4 = function (q, mu = 1, sigma = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(q < 0))
    stop(paste("y must be greater than 0 ", "\n", ""))
  sigma2 = exp(-sigma / mu)
  cdf = pweibull(q, shape = mu, scale = sigma2, lower.tail = lower.tail,
                 log.p = log.p)
  cdf
}


##########################################################
#Log-verossimilhança do modelo de mistura                                 #
##########################################################
Log_vero_mistura=function(vpar) {
  ph0 = vpar[1]
  ph1 = vpar[2]
  beta = vpar[-c(1:2)]
  linear=c(X%*%beta)
  theta=plogis(linear) # ligaÃ§Ã£o logistica 
  vf = dWEI4(y,mu=ph0,sigma=ph1)
  vF = pWEI4(y,mu=ph0,sigma=ph1) 
  Spop <- theta+(1-theta)*(1-vF)     
  fpop <- (1-theta)*vf 
  loglik  <- sum(status*log(fpop)+(1-status)*log(Spop))      
  loglik
}
##########################################################
#Log-verossimilhança do modelo de tempo de promoção                               #
##########################################################
Log_vero_p=function(vpar) {
  ph0 = (vpar[1])
  ph1 = vpar[2]
  beta = vpar[-c(1:2)]
  linear = c(X%*%beta)
  theta = exp(linear) 
  vf = dWEI4(y,mu=ph0,sigma=ph1)
  vF = pWEI4(y,mu=ph0,sigma=ph1) 
  Spop <- exp(-theta*(vF))     
  fpop <- theta*Spop*vf 
  loglik  <- sum(status*log(fpop)+(1-status)*log(Spop))      
  loglik
}
###############################################################
#Log-verossimilhança do modelo tempo de promoção com dispersão                                #
###############################################################
#########################BN############
Log_vero_b=function(vpar) {
  gamma=vpar[1]
  ph0 = vpar[2]
  ph1 = vpar[3]
  beta = vpar[-c(1:3)]
  linear=(X%*%beta)
  theta=exp(linear) # ligação log
  vf = dWEI4(y,mu=ph0,sigma=ph1)
  vF = pWEI4(y,mu=ph0,sigma=ph1) 
  Spop=(1+theta*gamma*vF)^(-1/gamma)
  fpop=theta*vf*(1+theta*gamma*vF)^(-1/gamma-1)
  loglik  <- sum(status*log(fpop)+(1-status)*log(Spop))      
  loglik
}
####
library(survival)

# Session -> Set Working Directory -> To Source File Locationl
dados <- read.table("melanoma.txt", header=TRUE, sep="\t")
head(dados)
y=dados$t   #years
status=dados$d1
x1 = dados$x1	 #treatment 0:observation and 1:interferon  
x2 = ifelse(dados$x2>48,1,0)
x3 = as.factor(dados$x3)  #nodule (levels:1,2,3,4)
x4 = as.factor(dados$x4) #sex (0:male and 1:female)
x5 = as.factor(dados$x5) 
x6 = dados$x6

#############################################
# Ajuste do modelo de mistura Weibull
##########################################
X=model.matrix(~1+x1+x2+x3+x4+x5+x6)
nc=ncol(X)

vpar1 = c(1,1,rep(0.1,nc))
Log_vero_mistura(vpar1)
fitm=optim(vpar1,Log_vero_mistura,control=list(fnscale=-1),method="L-BFGS-B",hessian = T,
           lower=c(0.0001,-Inf,rep(-Inf,(nc))),upper=c(rep(Inf,length(vpar1))))
EMVm=fitm$par
covm=solve(-fitm$hessian)
EPm=sqrt(diag(covm))

Estm=EMVm/EPm
p_valorm=2*(1-pnorm(abs(Estm)))
Saidam=cbind(EMVm, EPm,Estm,p_valorm)
rownames(Saidam)=c("phi_1","phi_2",paste("beta_", 0:(nc-1), sep = ""))
print(Saidam, digits = 3)
AIC_m=-2*fitm$value+2*length(EMVm)
library(xtable)
xtable(Saidam,digits=3)
#############################################
# Ajuste do modelo tempo de promoção  Weibull
##########################################

vpar1 = c(1,1,rep(-0.1,nc))
Log_vero_p(vpar1)
fitp=optim(vpar1,Log_vero_p,control=list(fnscale=-1),method="L-BFGS-B",hessian = T,
          lower=c(0.00001,-Inf,rep(-Inf,(nc))),upper=c(rep(Inf,length(vpar1))))
EMVp=fitp$par
covp=solve(-fitp$hessian)
EPp=sqrt(diag(covp))

Estp=EMVp/EPp
p_valorp=2*(1-pnorm(abs(Estp)))
Saidap=cbind(EMVp, EPp,p_valorp)
rownames(Saidap)=c("phi_1","phi_2",paste("beta_", 0:(nc-1), sep = ""))
print(Saidap, digits = 3)
AIC_p=-2*fitp$value+2*length(EMVp)
library(xtable)
xtable(Saidap,digits=3)
###################################################################
# Ajuste do modelo tempo de promoção  Weibull com Dispersão
##################################################################

vpar1 = c(0.5,1,-log(mean(y)),rep(1,nc))
fitb=optim(vpar1,Log_vero_b,control=list(fnscale=-1),method="L-BFGS-B",hessian = T,
           lower=c(0.0001,0.0001,-Inf,rep(-Inf,(nc))),upper=c(rep(Inf,length(vpar1))))

EMVb=fitb$par
covb=solve(-fitb$hessian)
EPb=sqrt(diag(covb))

Estb=EMVb/EPb
p_valorb=2*(1-pnorm(abs(Estb)))
Saidab=cbind(EMVb, EPb,p_valorb)
rownames(Saidab)=c("eta","phi_1","phi_2",paste("beta_", 0:(nc-1), sep = ""))
print(Saidab, digits = 3)
AIC_b=-2*fitb$value+2*length(EMVb) 
library(xtable)
xtable(Saidab,digits=3)
AIC=rbind(AIC_m,AIC_p,AIC_b)
colnames(AIC)=c("Críterio AIC")
rownames(AIC)=c(" Mistura","Tempo de promoção","Tempo de promoção-Dispersão")
xtable(AIC,digits = 3)
### ############################################################
#Ajuste do modelo tempo de promoção com dispersão reduzido
#############################################################
X=model.matrix(~1+x3)
nc=ncol(X)
vpar1 = c(0.5,1,-log(mean(y)),rep(1,nc))
fitbr=optim(vpar1,Log_vero_b,control=list(fnscale=-1),method="L-BFGS-B",hessian = T,
           lower=c(0.0001,0.0001,-Inf,rep(-Inf,(nc))),upper=c(rep(Inf,length(vpar1))))

EMVbr=fitbr$par
covbr=solve(-fitbr$hessian)
EPbr=sqrt(diag(covbr))

Estbr=EMVbr/EPbr
p_valorbr=2*(1-pnorm(abs(Estbr)))
Saidabr=cbind(EMVbr, EPbr,p_valorbr)
rownames(Saidabr)=c("eta","phi_1","phi_2",paste("beta_", 0:(nc-1), sep = ""))
print(Saidabr, digits = 3)
fitbr$value

xtable(Saidabr,digits=3)

TRV=-2*(fitbr$value-fitb$value)
pvalor=1-pchisq(TRV,5)
print(cbind(TRV,pvalor))
####################################
## Função a proporção de curados
###################################
p0=function(vpar,X) {
  gamma=vpar[1]
  ph0 = vpar[2]
  ph1 = vpar[3]
  beta = vpar[-c(1:3)]
  linear=c(X%*%beta)
  theta=exp(linear) # ligaÃ§Ã£o log
  p0=(1+theta*gamma)^(-1/gamma)
  p0
}
####################################
#EMV da p0 
####################################
X1=c(1,0,0,0)
X2=c(1,1,0,0)
X3=c(1,0,1,0)
X4=c(1,0,0,1)
XX=rbind(X1,X2,X3,X4)
pp0=p0(EMVbr,XX)
EP1=numeric()
for(j in 1:4){
g1=rbind(grad(p0,EMVbr,X=XX[j,]))
EP1[j]=sqrt(g1%*%covbr%*%t(g1))
}
LI=pp0-qnorm(0.975)*EP1
LS=pp0+qnorm(0.975)*EP1
saida=cbind(XX,pp0,LI,LS)
xtable(saida,digits = 3)


