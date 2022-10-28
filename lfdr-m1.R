##Implementation of 2d rejection region with a marginal local fdr
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/EM2v2.R")
library(ggplot2)
##Generate data:
size = 0.05
J = 1000
mu = 5
theta = -5
pi = c(0.7, 0.1, 0.1, 0.1)
var1 = 1
var2 = 2

obj = generate(J, pi, var1, var2, mu, theta)
samp = obj$sample
gamma = obj$gamma
tp = (gamma == 4)
tn = (gamma != 4)


##EM algorithm to estimate the parameters.
pi.init = c(0.8, 0.05, 0.05, 0.1)
parm = maximization(samp, pi.init)$par

##Will use EM algorithm to estimate the parameters later on, 
#for now, using the actual parameter values to compute the oracle rule.


lfdr2 = function(samp, parm)
{
  pi = parm[1:4]
  mu = parm[5]
  theta = parm[6]
  var1 = parm[7]
  var2 = parm[8]
  a = samp[,1]
  b = samp[,2]
  lfdra = vector()
  lfdrb = vector()
  for(i in 1:J)
  {
    t1 = (pi[1] + pi[3])*dnorm(a[i], mean = 0, sd = sqrt(var1))
    t2 = (pi[2] + pi[4])*dnorm(a[i], mean = mu, sd = sqrt(var1 + 1))
    lfdra[i] = t1/(t1 + t2)
    
    t3 = (pi[1] + pi[2])*dnorm(b[i], mean = 0, sd = sqrt(var2))
    t4 = (pi[3] + pi[4])*dnorm(b[i], mean = theta, sd = sqrt(var2 + 1))
    lfdrb[i] = t3/(t3 + t4)
  }
  return(list(lfdra = lfdra, lfdrb = lfdrb))
}

obj = lfdr2(samp, parm = parm)
lfdra = obj[[1]]
lfdrb = obj[[2]]


##Making a ggplot
hyp = if_else(tp, "non_null", "null")
data = data.frame(lfdra, lfdrb, hyp)
ggplot(data, aes(x = lfdra, y = lfdrb, color = hyp)) + geom_point()


#Function for computing mFDR
Q = function(i,j, lfdra, lfdrb, parm)
{
  pi = parm[1:4]
  st.lfdra = sort(lfdra)
  st.lfdrb = sort(lfdrb)
  Am = sum(st.lfdra[1:i])
  Bm = sum(st.lfdrb[1:j])
  Am.til = sum((1-st.lfdra)[1:i])
  Bm.til = sum((1-st.lfdrb)[1:j])
  if(i!=J && j!= J){
    Am.hat = sum((1-st.lfdra)[(i+1):J])
    Bm.hat = sum((1-st.lfdrb)[(j+1):J]) 
  }else{
    Am.hat = 0
    Bm.hat = 0
  }
  
  t1 = pi[1]/((pi[1]+pi[2])*(pi[1] + pi[3]))*Am*Bm
  t2 = pi[3]/((pi[1] + pi[3])*(pi[3] + pi[4]))*Am*Bm.til
  t3 = pi[2]/((pi[1] + pi[2])*(pi[2] + pi[4]))*Am.til*Bm
  t4 = pi[4]/((pi[2] + pi[4])*(pi[3] + pi[4]))*Am.hat*Bm.hat
  den = max(1,sum((lfdra < st.lfdra[i])*(lfdrb < st.lfdrb[j])))
  
  
  return(list(mfdr = (t1 + t2 + t3)/den, mfnr= t4/den))
}



#Finding cutoffs
cutoff = function(lfdra, lfdrb, q)
{
  Qmfdr = matrix(nrow = J, ncol = J)
  Qmfnr = matrix(nrow = J, ncol = J)
  G = matrix(nrow = J, ncol = J)
  for(i in 1:J)
  {
    for(j in 1:J)
    {
      obj= Q(i,j, lfdra, lfdrb, parm = parm)
      Qmfdr[i,j] = obj$mfdr
      Qmfnr[i,j] = obj$mfnr
      
      
    }
    print(i)
  }
  G = (Qmfdr <= size) + 0
if(sum(G) == 0) {
  rej = rep(0,J)
  print("Higher cutoff necessary")
}else{
  temp = (max(Qmfnr) -Qmfnr)*G
  ###Fix this tomorrow. You're looking for the minimum value among the 
  ###non zero values of temp. Simple, find a monotone decreasing transform of the mfdr
  cutoffs = which(temp == max(temp), arr.ind = TRUE)[1,]
  lambda = sort(lfdra)[cutoffs[1]]
  psi = sort(lfdrb)[cutoffs[2]]
  reject = (lfdra<=lambda)*(lfdrb <=psi)
}
}

reject_new = (lfdra<= 0.75)*(lfdrb <= 0.54)
sum(reject_new)
mfdr = sum(tn*reject_new)/sum(reject_new)
pow = sum(reject_new*tp)/sum(tp)
mfdr
pow
