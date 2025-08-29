# Limpeza do ambiente de trabalho
rm(list=ls(all=TRUE))

#Bibliotecas 
library(survminer)
library(survival)
library(ggfortify)


# Tempos e censuras fornecidas pelo exercicio:
temposCompath <- c(8,11,19,24,28,33,36,38,44,96,124,130,250,250,250)
censuraCompath <- c(1,1,1,0,1,1,0,1,1,1,1,1,1,0,0)

temposZena <- c(7,8,10,12,13,14,19,23,25,26,27,31,31,49,59,64,87,89,107,117,119,
                130,148,153,156,159,191,222,rep(250,16))
censuraZena <- c(rep(1,5),0,1,1,0,1,1,1,0,1,0,0,rep(1,12),rep(0,16))


# Vetor contendo os tempos para as duas drogas (concatenado)
tempos <- c(temposCompath,temposZena)
# Vetor contendo as censuras para as duas drogas (concatenado)
censura <- c(censuraCompath,censuraZena)

# Divisor para os grupos de estudo (grupo 1: Compath, grupo 2: Zena)
grupo <- c(rep(1,length(temposCompath)), rep(2,length(temposZena)))

#Teste Logrank:
fit <- survdiff(Surv(tempos,censura)~grupo)
fit

            #######SAIDA OBSERVADA#######


#         N Observed Expected (O-E)^2/E (O-E)^2/V
#grupo=1 15       11     6.76     2.656      3.38
#grupo=2 44       23    27.24     0.659      3.38

#   Chisq= 3.4  on 1 degrees of freedom, p= 0.07 


rm(list=ls(all=TRUE))

library(survminer)
library(survival)
library(ggfortify)

require(survival)
# 1, 1, 2, 4, 4, 6, 6, 7, 8, 9
# 9, 10, 12, 13, 14, 18, 19, 24, 26, 29
# 31+, 42, 45+ ,50+, 57, 60, 71+, 85+ 91
tempos<- c(1,1,2,4,4,6,6,7,8,9,9,10,12,13,14,18,19,24,26,29,31,42,45,
           50,57,60,71,85,91)
cens<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,1,0)
grupo = rep(1,length(tempos))
ekm<- survfit(Surv(tempos,cens)~1)
autoplot(ekm)