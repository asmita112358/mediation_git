##Implementation of 1d rejection region with a joint local fdr
library(combinat)
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/EM2v2.R")

##Function for computing cutoffs
lfdr = function(samp, parm)
{
  J = nrow(samp)
  a = samp[,1]
  b = samp[,2]
  pi = parm[1:4]
  mu = parm[5]
  theta = parm[6]
  var1 = parm[7]
  var2 = parm[8]
  a = samp[,1]
  b = samp[,2]
  
  lfdr = vector()
  for(i in 1:J)
  {
    t1 = pi[1]*dnorm(a[i],0,sqrt(var1))*dnorm(b[i],0,sqrt(var2)) + 
      pi[2]*dnorm(a[i], mu,sqrt(var1 + 1))*dnorm(b[i],0,sqrt(var2)) + 
      pi[3]*dnorm(a[i],0,sqrt(var1))*dnorm(b[i],theta,sqrt(var2 + 1))
    t2 = pi[4]*dnorm(a[i], mu,sqrt(var1 + 1))*dnorm(b[i],theta,sqrt(var2 + 1))
    lfdr[i] = t1/(t1 + t2)
    #rm(list = c(t1,t2))
  }
  st.lfdr<-sort(lfdr)
  k=1
  q = 0.05
  while(k<J && ((1/k)*sum(st.lfdr[1:k])) <= q){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<- lfdr<=lfdrk
  accept<- lfdr>lfdrk
  return(list(cutoff = lfdrk, rej = reject, k = k))
}





sim.fcn = function(J, pi, var1, var2, mu, theta, pi.init = c(0.8, 0.05, 0.05, 0.1))
{
  obj = generate(J, pi, var1, var2, mu, theta)
  samp = obj$sample
  gamma = obj$gamma
  tp = (gamma == 4)
  tn = (gamma != 4)
 
  ##Running EM algorithm to estimate the parameters
  
  parm = maximization(samp, pi.init)$par
  
  
  ##Checking which rearrangement of pi suits
  
  pi_est = parm[1:4]
  
  err = vector()
  object = permn(pi_est)
  for(i in 1:24)
  {
    err[i] = sum(abs((object[[i]]-pi)))
  }
  pi_new = object[[which.min(err)]]
  
  ##computing the joint local fdr
  
  ##samp is a JX2 matrix with the first col as a, second col b
  ##Parm is the parameter vector
  
  parm = c(pi_new, parm[5:8])
  
  obj = lfdr(samp, parm = parm)
  reject = obj$rej
  pow = sum(reject*tp)/sum(tp)
  mfdr = sum(tn*reject)/sum(reject)
  
  return(c(mfdr, pow))
  
  
}


##Dense null:
##Generate data:
J = 1000
mu = 5
theta = 5
pi = c(0.7, 0.1, 0.1, 0.1)
var1 = 0.5
var2 = 0.7

m = matrix(nrow = 100, ncol = 2)
for(r in 1:100)
{
 m[r,] = sim.fcn(J, pi, var1, var2, mu, theta, pi.init = c(0.8, 0.05, 0.05, 0.1))
 print(r) 
}

u = colMeans(m)
v = apply(m,2, var)
v = sqrt(v)

##Dense alternative
J = 1000
mu = 5
theta = 5
pi = c(0.3, 0.2, 0.2, 0.3)
var1 = 0.5
var2 = 0.7


m = matrix(nrow = 20, ncol = 2)
for(r in 1:20)
{
  m[r,] = sim.fcn(J, pi, var1, var2, mu, theta, pi.init = c(0.25, 0.25, 0.25, 0.25))
  print(r) 
}

u = colMeans(m)
v = apply(m,2, var)
