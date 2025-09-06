rm(list=ls(all=TRUE))

# Bibliotecas
library(survminer)
library(survival)
library(ggfortify)

##############################################################################

#a) Obtencao das curvas de sobrevivencia estimadas no N-A e por K-M

# Dados fornecidos pelo exercicio
tempos<- c(1,1,2,4,4,6,6,7,8,9,9,10,12,13,14,18,19,24,26,29,31,42,45,
           50,57,60,71,85,91)
censuras <- c(rep(1,20),0,1,0,0,1,1,0,0,1)


# Funcao do pacote survival para obter o estimador K-M
fit1 <- survfit(Surv(tempos,censuras)~1)
# Funcao do pacote survival para obter o estimadorN-A
fit2 <- survfit(coxph(Surv(tempos,censuras) ~ 1))


# Alguns links interessantes para comandos em graficos de sobrevivencia:

# https://rpkgs.datanovia.com/survminer/
# https://github.com/kassambara/survminer/issues/195
# https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf
# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_surv.html


# Utilizacao do pacote survminer para a obtencao das curvas de sobrevivencia

fit <- list("Keplan-Meyer" = fit1, "Nelson-Aalen" = fit2)

ggsurvplot(fit, combine = TRUE,         # Combinar curvas no grafico     
           legend.labs =
             c("Keplan-Meyer",          # Legendas para as curvas
               "Nelson-Aalen"), 
           conf.int = T,                # Apresenta o IC[95%]   
           #conf.int.style = "step",    # Estilo grafico para o IC      
           censor = TRUE,               # Mostrar censuras      
           palette = "jco",
           xlab = "Tempo t (semanas)",  # Nomes para os eixos
           ylab = "P(T>t)",
           ggtheme = theme_gray()       # Estilo de plotagem do ggplot
           #risk.table = TRUE,          # Appresentar tabela de risco
           #risk.table.col = "strata",  # Cores na tabela de risco
           #risk.table.height = 0.25    # Altura da tabela de risco
           )

##############################################################################


#b) Obtencao de uma estimativa para o TMV a partir da estimacao K-M:

# Pode-se obter os valores de t onde há degraus na funcao de sobrevivencia e 
# as respectivas probabilidades de sobrevivencia para cada patamar a partir
# das colunas da funcao summary.

tempos <- summary(fit1)[[2]] #segunda coluna apresenta os tempos 
sobrev <- summary(fit1)[[6]] #sexta coluna apresenta P(T>t) para cada patamar


length(tempos)
length(sobrev)

#Funcao para a obtencao do TMV:
tempo_medio_est <- function(S,t){ # recebe S e t
  
  tmv <- S[1] # Valor inicial para o processo iterativo (t1)
  
  for(i in 1:(length(t)-1)){ # Soma de produtos do segundo termo do TMV estimado
    tmv <- tmv + S[i] * (t[i+1] - t[i])
  }
  
  return(tmv) # retorna o valor obtido
}


tempo_medio_est(sobrev,tempos) #calculando-se para o problema 
# Resposta: 30.39655 semanas


##############################################################################


#c) Abaixo temos a rotina para a obtencao do erro padrao 

#O valor r corresponde de falhas (observacoes nao censuradas)
r = length(censuras[censuras == 1]) 

# Obtencao dos dados da curva de sobrevivencia como no item anterior
tempos = summary(fit1)[[2]] 
sobrev = summary(fit1)[[6]]

# Utilizando-se dos metodos do objeto fit1, pode-se obter as seguintes listas:
n_ur = fit1$n.risk #lista de unidades em risco ate determinado t
n_uf = fit1$n.event#lista de unidades que falharam ate determinado t

# A correspondencia entre esses valores e os tempos podem ser verificadas na
# tabela gerada pelo comando summary(fit1)
  
var_tmv = 0#inicializacao do valor da variancia estimada do tempo medio de vida
  
A = rep(0, length(sobrev)) # inicializacao do vetor de valores de Aj
  
#Laco para a obtencao dos valores dos termos Aj 
for(i in 1:(length(sobrev) - 1) ){
  A[i] = sobrev[i]*(tempos[i + 1] - tempos[i])
}
  
#Laco para a obtencao da somatoria de Aj^2/nj(nj-dj), de j = 1 ate j = r-1
for(i in 1:(length(sobrev) - 1)){
  var_tmv = 
    var_tmv + sum(A[i:length(A2)])^2/(n_ur[i] * (n_ur[i] - n_uf[i]))
  }
  
var_tmv = (r/(r-1))*var_tmv #obtencao do valor final para a variancia do 
# estimador do tmv

  
ep_tmv <- sqrt(var_tmv) #obtencao do erro padrao estimado do estimador de tmv

ep_tmv

# Resposta: 5.373217










fit2<-survfit(Surv(tempos, censuras)~1, se.fit = TRUE)
fit3=summary(fit2)

cj=c(rep(0,15),1,2,0,2,0)
atuarial = data.frame(tempo=fit3$time, EmRisco=fit3$n.risk, 
                      Eventos=fit3$n.event, Censura=cj)
atuarial['qj']=
  1-atuarial['Eventos']/(atuarial['EmRisco']-0.5*atuarial['Censura'])

atuarial$sobrevivencia=atuarial$qj

for (i in 2:20) {
  atuarial$sobrevivencia[i] = 
    atuarial$sobrevivencia[i]*atuarial$sobrevivencia[i-1]  
} 

print(atuarial)
