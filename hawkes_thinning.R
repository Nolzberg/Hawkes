## Hawkes methode 1

# Packages ####

#install.packages("latex2exp")
library(latex2exp)

# Fonctions #### 

phi = function(x){
  alpha = 1.5
  delta = 2
  return(alpha*exp(-delta*x))
}

intensite = function(t){
  l = lambda0
  Tk = Ti[Ti < t]
  if (length(Tk) != 0){
    for (i in 1:length(Tk)){
      x = t - Tk[i]
      l = l + phi(x)
    }
  }
  return(l)
}

intensite_2para = function(t,Ti){
  l = lambda0
  if (length(Ti) != 0){
    for (i in 1:length(Ti)){
      if (Ti[i] != 0){
        x = t - Ti[i]
        l = l + phi(x)
      }
    }
  }
  return(l)
}

Hawkes_th = function(N){
  Tn = rep(0,N)
  t = 0
  i = 1
  while(i <= N){
    lambda_etoile = intensite_2para(t,Tn)
    u = runif(1)
    tau = tau = -log(u)/lambda_etoile
    t = t + tau
    s = runif(1)
    lambda = intensite_2para(t,Tn)
    if(s <= lambda/lambda_etoile){
      Tn[i] = t
      i = i + 1 
    }
  }
  return(c(0,Tn))
}

# Fonctions d'affichage ####

plot_intensite = function(t,Tn){
  tim = seq(0,t,by = 0.001)
  val = rep(NA,length(tim))
  d = 1
  for (i in tim){
    Ti = Tn[Tn < i]
    val[d] = intensite_2para(i,Ti)
    d = d + 1
  }
  plot(tim,val,type = 'l',main = TeX("Intensité du processus de Hawkes"), xlab = "t",ylab = TeX("$\\lambda (t)$"))
}

plot_Nt = function(t,Tn){
  tim = seq(0,t,by = 0.001)
  val = rep(NA,length(tim))
  d = 1
  for (i in tim){
    Ti = Tn[Tn < i & Tn != 0]
    val[d] = length(Ti)
    d = d + 1
  }
  plot(tim,val,type = 'l',main = TeX("Simulation d'un processus de Hawkes (méthode thinning)"), xlab = "t",ylab = TeX("$N(t)$"))
}

# Tests ####

N = 10
l = N + 1
lambda0 = 0.9

Ti = Hawkes_th(N)

par(mfrow=c(2,1))
plot_Nt(Ti[N+1]+0.3,Ti)
plot_intensite(Ti[l]+0.3,Ti)

