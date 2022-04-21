###
# Hawkes multivarié pour m neurones
###

## Packages ####

#install.packages("latex2exp")
library(latex2exp)

## Fonctions ####

Intensite = function(m,t,Nk,Y,lambda0,a,delta){
  l = a[m] + (lambda0[m] - a[m])*exp(-delta[m]*t)
  n = length(Nk)
  for(k in 1:n){
    Nk[[k]] = Nk[[k]][Nk[[k]] < t]
    if (length(Nk[[k]]) > 1){
      for (i in 2:length(Nk[[k]])){
        x = t - Nk[[k]][i]
        l = l + Y[[(m-1)*n+k]][i]*exp(-delta[m]*x)
      }
    }
  }
  return(l)
}

simu_multi = function(lambda0,a,delta,beta,N){
  
  Tn = rep(NA,N)
  Tn[1] = 0
  
  m = length(a)
  nb = rep(0,m*m)
  n_k = rep(0,m)
  Nk = list()
  Nk = append(Nk,n_k)
  Y = list()
  Y = append(Y,nb)
  
  lambda = matrix(c(NA,lambda0[1]),ncol = 2)
  
  for(i in 2:m){
    lambda = rbind(lambda,c(NA,lambda0[i]))
  }
  
  for(k in 1:N){
    s = c(NA,NA)
    for (i in 1:m){
      s0 = -log(runif(1))/a[i]
      u1 = runif(1)
      d = 1 + delta[i]*log(u1)/(lambda[i,2] - a[i])
      if(d > 0){
        s1 = -log(d)/delta[i]
        s[i] = min(s0,s1)
      }else{s[i] = s0}
    }
    
    w = min(s)
    l = which.min(s)
    Tn[k+1] = Tn[k] + w
    
    for(o in 1:m){
      if(l == o){
        Nk[[o]] = append(Nk[[o]],Tn[k+1])
      }
    }
    
    for(j in 1:m){
      lambda[j,1] = (lambda[j,2] - a[j])*exp(-delta[j]*(Tn[k+1]-Tn[k])) + a[j]
      y = rexp(1,beta[j,l])
      lambda[j,2] = lambda[j,1] + y
      Y[[(j-1)*m+l]] = append(Y[[(j-1)*m+l]],y)
    }
  }
  return(list('Tn' = Tn,'T_Nk' = Nk,'Yk' = Y,'Nombre de neurones' = m))
}

plot_Nt = function(Nk){
  m = length(Nk)
  col = rainbow(m)
  ma = 0
  for(i in 1:m){
    maxi = max(length(Nk[[i]]))
    if(maxi >= ma){
      ma = maxi
    }
  }
  plot(Nk[[1]],0:(length(Nk[[1]])-1),type='s',col = col[1], ylim = c(0,ma),main = TeX("Processus de Hawkes avec n neurones"), xlab = "t",ylab = TeX("$\\N_t$"))
  for(p in 2:m){
    lines(Nk[[p]],0:(length(Nk[[p]])-1),type='s',col = col[p])
  }
}

plot_inte = function(Tn,Nk,Y,lambda0,a,delta){
  m = length(Nk)
  time = seq(0,Tn[N+1],length.out = 500)
  l = length(time)
  val = matrix(rep(NA,m*l),nrow = m)
  d = 1
  for(p in time){
    for(k in 1:m){
      val[k,d] = Intensite(k,p,Nk,Y,lambda0,a,delta)
    }
    d = d + 1
  }
  col = rainbow(m)
  ma = max(val)
  plot(time,val[1,],type='l',col = col[1], ylim = c(0,ma),main = TeX("Intensité du processus de Hawkes avec n neurones"), xlab = "t",ylab = TeX("$\\lambda (t)$"))
  for(p in 2:m){
    lines(time,val[p,],type='l',col = col[p])
  }
  return(list('Temps' = time,'Intensite' = val))
}

plot_champMoy = function(temps,val){
  valmoy = colMeans(val)
  valmax = apply(val,2,max)
  valmin = apply(val,2,min)
  plot(temps,valmoy,type = 'l',col = 'black',ylim = c(0,max(val)),main = TeX("Champ Moyen de l'intensité du processus de Hawkes avec 10 neurones"), xlab = "t",ylab = TeX("$\\lambda (t)$"))
  lines(temps,valmax,type = 'l',col = 'red')
  lines(temps,valmin,type = 'l',col = 'lightblue')
  legend('topleft',legend = c('Champ Moyen','Max','Min'), col = c('black','red','lightblue'),lty = c(1,1,1),cex = 0.5)
}

## Test 3 neurones ####

lambda0 = c(0.7,0.7,0.7)
a = c(0.4,0.3,0.3)
delta = c(0.8,1.0,1.0)
beta = matrix(c(2.5,4.0,8.0,8.0,2.0,4.0,16.0,8.0,2.5),nrow = 3,byrow = TRUE)
N = 2000

test_3neurones = simu_multi(lambda0,a,delta,beta,N)

m3 = test_3neurones$`Nombre de neurones`
Y3 = test_3neurones$Yk
Nk3 = test_3neurones$T_Nk
Tn3 = test_3neurones$Tn

par(mfrow = c(1,1))
plot(Tn3,0:(length(Tn3)-1),type = 's')

par(mfrow = c(2,1))
plot_Nt(Nk3)
plot_inte(Tn3,Nk3,Y3,lambda0,a,delta)

## Test 10 neurones ####

lambda010 = rep(0.7,10)
a10 = c(0.4,0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3,0.4)
delta10 = c(0.8,1.0,1.0,0.8,1.0,0.8,0.8,1.0,0.8,1.0)
beta10 = matrix(sample(10:20,size = 100,replace = TRUE),nrow = 10,byrow = TRUE)
N = 3000

test_10neurones = simu_multi(lambda010,a10,delta10,beta10,N)
m10 = test_10neurones$`Nombre de neurones`
Y10 = test_10neurones$Yk
Nk10 = test_10neurones$T_Nk
Tn10 = test_10neurones$Tn

par(mfrow = c(2,1))
test = plot_inte(Tn10,Nk10,Y10,lambda010,a10,delta10)
plot_champMoy(test$Temps,test$Intensite)

par(mfrow = c(1,1))
plot(Tn10,0:(length(Tn10)-1),type = 's')

## Test 100 neurones ####

nb = 1000

lambda0100 = rep(0.7,nb)
a100 = sample(3:4,size = nb,replace = TRUE)/10
delta100 = sample(8:10,size = nb,replace = TRUE)/10
beta100 = matrix(sample(100:500,size = nb^2,replace = TRUE),nrow = nb,byrow = TRUE)
for(i in 1:nb){
  beta100[i,i] = sample(3:5,size = 1)
}
N = 10000

test_100neurones = simu_multi(lambda0100,a100,delta100,beta100,N)
m100 = test_100neurones$`Nombre de neurones`
Y100 = test_100neurones$Yk
Nk100 = test_100neurones$T_Nk
Tn100 = test_100neurones$Tn

par(mfrow = c(1,1))
plot(Tn100,0:(length(Tn100)-1),type = 's')

par(mfrow = c(3,1))
plot_Nt(Nk100)
test100 = plot_inte(Tn100,Nk100,Y100,lambda0100,a100,delta100)
plot_champMoy(test100$Temps,test100$Intensite)

## Test 2 neurones ####

lambda0_2 = c(0.7,0.7)
a_2 = c(0.4,0.6)
delta_2 = c(0.8,1.0)
beta_2 = matrix(c(1.5,4.0,8.0,2.0),nrow = 2,byrow = TRUE)

N = 5000

test_2neurones = simu_multi(lambda0_2,a_2,delta_2,beta_2,N)
m2 = test_2neurones$`Nombre de neurones`
Y2 = test_2neurones$Yk
Nk2 = test_2neurones$T_Nk
Tn2 = test_2neurones$Tn

par(mar = c(2,2),mfrow = c(2,1))
plot_Nt(Nk2)
test2 = plot_inte(Tn2,Nk2,Y2,lambda0_2,a_2,delta_2)

par(mfrow = c(1,1))
plot(Tn2,0:(length(Tn2)-1),type = 's')

