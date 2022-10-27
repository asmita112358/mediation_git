##Implementation of 1d rejection region with a joint local fdr
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/EM2v2.R")
##Generate data:
J = 1000
mu = 3
theta = -3
pi = c(0.7, 0.1, 0.1, 0.1)
var1 = 1
var2 = 2

obj = generate(J, pi, var1, var2, mu, theta)
samp = obj$sample
gamma = obj$gamma
tp = (gamma == 4)
tn = (gamma != 4)
##Running EM algorithm to estimate the parameters
pi.init = c(0.8, 0.05, 0.05, 0.1)
parm = maximization(samp, pi.init)$par

##computing the joint local fdr

##samp is a JX2 matrix with the first col as a, second col b
##parm = the parameter in the order it's returned from the maximization function, barring the pi's.

##Check this every time, or come up with a method to check which is closest

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

obj = lfdr(samp, parm = c(pi, mu, theta, var1, var2))
reject = obj$rej
pow = sum(reject*tp)/sum(tp)
mfdr = sum(tn*reject)/sum(reject)

pow
mfdr
