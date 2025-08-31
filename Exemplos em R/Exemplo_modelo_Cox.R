laringe<-read.table("https://docs.ufpr.br/~giolo/Livro/ApendiceA/laringe.txt", h=T)     
attach(laringe)
library(survival)
fit4<-coxph(Surv(tempos,cens) ~ factor(estagio) + idade + factor(estagio)*idade, data=laringe, 
            x = T, method="breslow")
#######################################
# An??lise residual                    #
#######################################
resid(fit4, type="scaledsch")
cox.zph(fit4,transform = "identity")
par(mfrow=c(1,3))
plot(cox.zph(fit4))
#######################################
#Figuras da Sobreviv??ncia e Risco     #
#######################################

Ht<-basehaz(fit4, centered=F)
tempos<-Ht$time
H0<-Ht$hazard
S0<- exp(-H0)
round(cbind(tempos,S0,H0), digits=5)

tt<-sort(tempos)
aux1<-as.matrix(tt)
n<-nrow(aux1)
aux2<-as.matrix(cbind(tempos,S0))
S00<-rep(max(aux2[,2]),n)
for(i in 1:n){
  if(tt[i]> min(aux2[,1])){
    i1<- aux2[,1]<= tt[i]
    S00[i]<-min(aux2[i1,2])}}
ts0<-cbind(tt,S00)
ts0

b<-fit4$coefficients
id<-50
st1<-S00^(exp(b[4]*id))                # S(t|x) estagio I   e idade = 50 anos
st2<-S00^(exp(b[1]+((b[4]+b[5])*id)))  # S(t|x) estagio II  e idade = 50 anos
st3<-S00^(exp(b[2]+((b[4]+b[6])*id)))  # S(t|x) estagio III e idade = 50 anos
st4<-S00^(exp(b[3]+((b[4]+b[7])*id)))  # S(t|x) estagio IV  e idade = 50 anos
id<-65
st11<-S00^(exp(b[4]*id))               # S(t|x) estagio I   e idade = 65 anos
st21<-S00^(exp(b[1]+((b[4]+b[5])*id))) # S(t|x) estagio II  e idade = 65 anos
st31<-S00^(exp(b[2]+((b[4]+b[6])*id))) # S(t|x) estagio III e idade = 65 anos
st41<-S00^(exp(b[3]+((b[4]+b[7])*id))) # S(t|x) estagio IV  e idade = 65 anos

par(mfrow=c(1,2))
plot(tt,st1, type="s", ylim=range(c(0,1)), xlab="Tempos", ylab="S(t|x)", lty=1)
lines(tt,st2, type="s", lty=2)
lines(tt,st3, type="s", lty=3)
lines(tt,st4, type="s", lty=4)
legend(0, 0.2, lty=c(1,2,3,4), c("estágio I","estágio II","estágio III","est??gio IV"), 
       lwd=1, bty="n", cex=0.7)
title("Idade = 50 anos")

plot(tt,st11, type="s", ylim=range(c(0,1)), xlab="Tempos", ylab="S(t|x)", lty=1)
lines(tt,st21,type="s",lty=2)
lines(tt,st31,type="s",lty=3)
lines(tt,st41,type="s",lty=4)
legend(0, 0.2, lty=c(1,2,3,4), c("estagio I","estagio II","estagio III","estagio IV"), 
       lwd=1, bty="n", cex=0.7)
title("Idade = 65 anos")

###
Ht1<- -log(st1)
Ht2<- -log(st2)
Ht3<- -log(st3)
Ht4<- -log(st4)
Ht11<- -log(st11)
Ht21<- -log(st21)
Ht31<- -log(st31)
Ht41<- -log(st41)
par(mfrow=c(1,2))
plot(tt, Ht1, type="s", ylim=range(c(0,4)), xlab="Tempos", ylab="Risco Acumulado", lty=1)
lines(tt, Ht2, type="s", lty=2)
lines(tt, Ht3, type="s", lty=3)
lines(tt, Ht4, type="s", lty=4)
legend(0.5,3.5, lty=c(1,2,3,4), c("est??gio I","est??gio II","est??gio III","est??gio IV"),
       lwd=1, bty="n", cex=0.7)
title("Idade = 50 anos")
plot(tt, Ht11, type="s", ylim=range(c(0,4)), xlab="Tempos", ylab="Risco Acumulado", lty=1)
lines(tt, Ht21, type="s", lty=2)
lines(tt, Ht31, type="s", lty=3)
lines(tt, Ht41, type="s", lty=4)
legend(0.5,3.5, lty=c(1,2,3,4), c("estadio I","est??gio II","est??gio III","est??gio IV"),
       lwd=1, bty="n", cex=0.7)
title("Idade = 65 anos")
