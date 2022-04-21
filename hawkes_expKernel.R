###
# TER : Simulation processus de Hawkes
###

## Packages ####

#install.packages("latex2exp")
library(latex2exp)

## Fonctions ####

simu_expKernel = function(N,a,delta,beta,lambda0,Yk){
  Tn = rep(NA,N)
  lambda = c(NA,lambda0)
  Tn[1] = 0
  for (i in 1:N){
    s0 = -log(runif(1))/a
    u1 = runif(1)
    d = 1 + delta*log(u1)/(lambda[2] - a)
    if(d > 0){
      s1 = -log(d)/delta
      tau = min(s0,s1)
    }else{tau = s0}
    Tn[i+1] = Tn[i] + tau
    lambda[1] = (lambda[2] - a)*exp(-delta*tau) + a
    lambda[2] = lambda[1] + Yk[i]
  }
  return(Tn)
}

g = function(x){lambda0 + 0.9*exp(-1.2*x)}

intensite = function(t){
  l = lambda0
  Tk = Tn[Tn < t]
  if (length(Tk) != 0){
    for (i in 1:length(Tk)){
      x = t - Tk[i]
      l = l + g(x)
    }
  }
  return(l)
}

## Variables et Tests ####

a = 0.9
delta = 1.0
beta = 1.2
lambda0 = 0.9
N = 50
Y = rexp(N,beta)

test2 = simu_expKernel(N,a,delta,beta,lambda0,Y)
plot(test2,0:N,type = 's')

