
rm(list=ls(all=TRUE))


t1=c(8,109,219,423,514,973,1013,1100,1246,1503,1607,rep(1680,29))
d1=c(rep(1,11),rep(0,29))
x1=rep(75,40)
t2=c(2,67,155,160,161,196,204,308,518,523,551,560,639,743,802,855,943,1236,1266,1287,1307,
     1351,1379,1321,1379,1521,1607,1639,1653,rep(1680,11))
d2=c(rep(1,27),rep(0,13))
x2=rep(95,40)

t3=c(1,11,11,12,24,32,59,63,75,81,84,86,130,133,135,142,142,195,202, 202,215,224,224,239,
     299,352,362,365,382,390,471,570,584,600,730,850,925,958,1413,1608)
d3=c(rep(1,39),rep(0,1))
x3=rep(115,40)
library(survival)
library(msm)

mKM_1=survfit(Surv(t1, d1) ~ 1, se.fit = FALSE)
mKM_2=survfit(Surv(t2, d2) ~ 1, se.fit = FALSE)
mKM_3=survfit(Surv(t3, d3) ~ 1, se.fit = FALSE)

par(mfrow=c(1,1))

### Avaliar Weibull
plot(log(mKM_1$time),log(-log(mKM_1$surv)) ,xlab="log (t)",ylab="log(-log(S(t)))", 
     ylim=c(-4,2),lty=1, lwd=2.5,pch=1,col=1)
points(log(mKM_2$time),log(-log(mKM_2$surv)) ,xlab="log (t)",ylab="log(-log(S(t)))", 
       lty=2, lwd=2.5,pch=2,col=2)
points(log(mKM_3$time),log(-log(mKM_3$surv)) ,xlab="log (t)",ylab="log(-log(S(t)))",lty=3,
       lwd=2.5,pch=3,col=3)

legend("topleft",c("75 ","95","115"), bty="n", pch=1:3,col=1:3,lty=1:3,lwd=2.5)

# LOGNormal
plot(log(mKM_1$time),qnorm(mKM_1$surv) ,xlab="log (t)",ylab=expression(Phi^-1 (S(t))), 
     ylim=c(-4,2),lty=1, lwd=2.5,pch=1,col=1)
points(log(mKM_2$time),qnorm(mKM_2$surv) ,xlab="log (t)",ylab=expression(Phi^-1 (S(t))), 
       lty=2, lwd=2.5,pch=2,col=2)
points(log(mKM_3$time),qnorm(mKM_3$surv) ,xlab="log (t)",ylab=expression(Phi^-1 (S(t))),lty=3,
       lwd=2.5,pch=3,col=3)
legend("bottomleft",c("75 ","95","115"), bty="n", pch=1:3,col=1:3,lty=1:3,lwd=2.5)

## LOGLOGISTICO
plot(log(mKM_1$time),log((1-mKM_1$surv)/mKM_1$surv), xlab="log(Tempo)",ylab="log((1-S(t))/S(t))", 
     ylim=c(-4,2),lty=1, lwd=2.5,pch=1,col=1)
points(log(mKM_2$time),log((1-mKM_2$surv)/mKM_2$surv), xlab="log(Tempo)",ylab="log((1-S(t))/S(t))", 
       lty=2, lwd=2.5,pch=2,col=2)
points(log(mKM_3$time),log((1-mKM_3$surv)/mKM_3$surv), xlab="log(Tempo)",ylab="log((1-S(t))/S(t))",
       lty=3,       lwd=2.5,pch=3,col=3)
legend("topleft",c("75 ","95","115"), bty="n", pch=1:3,col=1:3,lty=1:3,lwd=2.5)



t=c(t1,t2,t3)
delta=c(d1,d2,d3)
x=c(x1,x2,x3)
xx1=1000/(273.16+x)
mKM=survfit(Surv(t, delta) ~ x, se.fit = FALSE)
plot(mKM, xlab="Tempo (horas)",ylab="Fun??o de confiabilidade", lty=1:3, lwd=2.5,pch=1:3,col=1:3)
legend("bottomleft",c("75 ","95","115"), bty="n", pch=1:3,col=1:3,lty=1:3,lwd=2.5)

plot(x,t, xlab="Temperatura (oC)",ylab="Tempo de falha (hrs)")

fit.wei=survreg(Surv(t, delta)~xx1, dist="weibull",score=T)
summary(fit.wei)
fit.ln=survreg(Surv(t, delta)~xx1, dist="lognormal")
summary(fit.ln)
fit.ll=survreg(Surv(t, delta)~xx1, dist="loglogistic")
summary(fit.ll)
#dados=cbind(t,delta,x)
#write.table(dados,"d:\\marron.txt",row.names=F)
par(mfrow=c(1,1))
## An?lise de  Residuos-Weibull
SW=exp(-(t/exp(fit.wei$linear.predictors))^(1/fit.wei$scale))
eW=-log(SW)
eKMW=survfit(Surv(eW, delta) ~ 1, se.fit = FALSE)
plot(exp(-eKMW$time),eKMW$surv, xlab="S(e)-Exponencial",ylab="S(e)-KM",lwd=3, main="Weibull")
plot(eKMW,ylab="S(e)",xlab="Residuo-Cox-Snell",lty=1,lwd=2,main="Weibull")
curve(exp(-x),0,max(eKMW$time), add=T,lwd=2, lty=2)
legend("topright", c("KM","Exponencial"), bty="n", lty=1:2)
## An?lise de  Residuos-logNormal
SLN=pnorm(-(log(t)-fit.ln$linear.predictors)/fit.ln$scale)
eLN=-log(SLN)
eKML=survfit(Surv(eLN, delta) ~ 1, se.fit = FALSE)
plot(exp(-eKML$time),eKML$surv, xlab="S(e)-Exponencial",ylab="S(e)-KM",lwd=3,main="LogNormal")
plot(eKML,ylab="S(e)",xlab="Residuo-Cox-Snell",lty=1,lwd=2,main="LogNormal")
curve(exp(-x),0,max(eKML$time), add=T,lwd=2, lty=2)
legend("topright", c("KM","Exponencial"), bty="n", lty=1:2)
## An?lise de  Residuos-loglogistic
SLL=(1+(fit.ll$linear.predictors*t)^(1/fit.wei$scale))^(-1)
eLL=-log(SLN)
eKMLL=survfit(Surv(eLL, delta) ~ 1, se.fit = FALSE)
plot(exp(-eKMLL$time),eKMLL$surv, xlab="S(e)-Exponencial",ylab="S(e)-KM",lwd=3,main="LogLogistico")
plot(eKMLL,ylab="S(e)",xlab="Residuo-Cox-Snell",lty=1,lwd=2,main="LogLogistico")
curve(exp(-x),0,max(eKMLL$time), add=T,lwd=2, lty=2)
legend("topright", c("KM","Exponencial"), bty="n", lty=1:2)
###
betas=fit.wei$coefficients
alpha=1/fit.wei$scale
cov=fit.wei$var
ep=sqrt(diag(cov))
emv=c(betas,fit.wei$icoef[2])
## Delta method
dm=deltamethod(list(~1/exp(x3)), emv, cov)
mv=c(alpha,betas)
ep1=c(dm,ep[1:2])
z=mv/ep1
saida=cbind(mv,ep1,z)
colnames(saida)=c("Estimativa","Erro padr?o", "Estat.")
rownames(saida)=c("alpha","beta_0", "beta_1.")
(round(saida,4))
### Tempo medio
c=55
xc=1000/(273.16+c)
ET=exp(betas[1]+betas[2]*xc)*gamma((fit.wei$scale)+1)
EPE=deltamethod(list(~exp(x1+x2*xc)*gamma(exp(x3)+1)), emv, cov)
LIE=max(0,ET-1.64*EPE)
LSE= ET+ 1.64*EPE
ESE=cbind(ET,EPE,LIE,LSE)
colnames(ESE)=c("Valor Esperado", "Erro padr?o","LI", "LS")
ESE
## A media em fun??o do niveis de estresse
et=function(x){
  
  xc=1000/(273.16+x)
  ee= exp(betas[1]+betas[2]*xc)*gamma((fit.wei$scale)+1)
  ee
}
curve(et(x),45,120, xlab="Temperatura(oC)",ylab="Tempo m?dio (hrs)",lty=1)
segments(55,0,55,et(55),col="red",lty=2)
segments(0,et(55),55,et(55),col="red",lty=2)
### A mediana em fun??o dos n?veis de estresse
em=function(x){
  
  xc=1000/(273.16+x)
  ee= exp(betas[1]+betas[2]*xc)*log(2)^((fit.wei$scale))
  ee
}
curve(em(x),45,120, xlab="Temperatura(oC)",ylab="Tempo mediano")
segments(55,0,55,em(55),col="red",lty=2)
segments(0,em(55),55,em(55),col="red",lty=2)
MDE=deltamethod(list(~exp(x1+x2*xc)*log(2)^(exp(x3))), emv, cov)
md=em(55)
LIMd=max(0,md-1.64*MDE)
LSMd= md+ 1.64*MDE
ESMd=cbind(md,MDE,LIMd,LSMd)
colnames(ESMd)=c("Mediana", "Erro padr?o","LI", "LS")
ESMd
## Quantil fun??o
qwar=function(x,p=0.5){
  
  xc=1000/(273.16+x)
  ee= exp(betas[1]+betas[2]*xc)*(-log(1-p))^((fit.wei$scale))
  ee
}
curve(qwar(x,p=0.05),45,120, xlab="Temperatura(oC)",ylab="Tempo quantil 5%",col="black")
segments(55,0,55,qwar(55,0.05),col="red",lty=2)
segments(0,qwar(55,0.05),55,qwar(55,0.05),col="red",lty=2)

EpQ=deltamethod(list(~exp(x1+x2*xc)*(-log(1-0.05))^(exp(x3))), emv, cov)
q=qwar(55,0.05)
LIq=max(0,q-1.64*EpQ)
LSq= q+ 1.64*EpQ
Eq=cbind(q,EpQ,LIq,LSq)
colnames(Eq)=c("Quatil 5%", "Erro padr?o","LI", "LS")
Eq

## Reability function

Swar=function(x){
  
  xc=1000/(273.16+55)
  ee= exp(-(x/exp(betas[1]+betas[2]*xc))^(1/(fit.wei$scale)))
  ee
}

curve(Swar(x),0,150000, xlab="Tempo(hrs)",ylab="Fun??o de confiabilidade")
segments(1825,0,1825,Swar(1825),col="red",lty=2)
segments(0,Swar(1825),1825,Swar(1825),col="red",lty=2)
### Erro padr?o para S(to)
EPS=deltamethod(list(~exp(-(1825/exp(x1+x2*xc))^(1/exp(x3)))), emv, cov)
S_0=Swar(1825)
LIs=max(0,S_0-1.64*EPS)
LSs= S_0+ 1.64*EPS
Es=cbind(S_0,EPS,LIs,LSs)
colnames(Es)=c("S(t_0)", "Erro padr?o","LI", "LS")
Es

### Erro padr?o para 1-S(to), to=1825
EPF=deltamethod(list(~1-exp(-(1825/exp(x1+x2*xc))^(1/exp(x3)))), emv, cov)
F_0=1-Swar(1825)
LIF=max(0,F_0-1.64*EPF)
LSF= F_0+ 1.64*EPF
EF=cbind(F_0,EPF,LIF,LSF)
colnames(EF)=c("F(t_0)", "Erro padr?o","LI", "LS")
EF

Resultado=round(rbind(ESE/1825,ESMd/1825,Eq/1825,Es,EF),4)
rownames(Resultado)=c("Media","Mediana","Quantil 5%", "S(to)","F(t_0)")
colnames(Resultado)=c("Estimava","Erro padr?o", "LI","LS")
Resultado


