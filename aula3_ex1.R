rm(list=ls(all=TRUE))

Rs <- function(x){
  return(2*exp(-0.06*x)-exp(-0.09*x))
  }

R <- function(x){
  return(exp(-0.03*x))
  }

curve(Rs,0,250, col = 'black', ylab = "Rs(x) e R(x)")
curve(R, col = 'red', add = TRUE)

legend("topright",c("Rs(x)", "R(x)"), col =c("black", "red"),lty = 1, bty = "n")

# Outro exemplo
Rs <- function(x){
  return(2*exp(-0.2*x)-exp(-0.4*x))
}

R <- function(x){
  return(exp(-0.1*x))
}

curve(Rs,0,250, col = 'black', ylab = "Rs(x) e R(x)", xlim = c(0,75))
curve(R, col = 'firebrick', add = TRUE)
abline(v = 2,)

legend("topright",c("Rs(x)", "R(x)"), col =c("black", "firebrick"),
       lty = 1, bty = "n")


library(ggplot2)

p <- ggplot() + xlim(0,70)

p +
  geom_function(aes(colour = "R(t)"), fun = R, size = 1) +
  geom_function(aes(colour = "Rs(t)"), fun = Rs, size = 1) +
  labs(title = "Funções de Confiabilidade") +
  #scale_color_manual(name = "Functions",
   #                  values = c("deeppink4", "deepskyblue4"), 
    #                 labels = c("R(t)", "Rs(t)"))+
  xlab("t") +
  ylab("P(T>t)") + 
  
  theme_grey(base_size = 20)




# Outro exemplo

# Outro exemplo
Rs <- function(x){
  return(exp(-0.07*x))
}

R <- function(x){
  return(exp(-0.04*x^0.5))
}

curve(Rs,0,250, col = 'black', ylab = "Rs(x) e R(x)", xlim = c(0,1000))
curve(R, col = 'firebrick', add = TRUE)


legend("topright",c("Rs(x)", "R(x)"), col =c("black", "firebrick"),
       lty = 1, bty = "n")
