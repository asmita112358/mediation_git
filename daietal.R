
library(DACT)
library(HDMT)
null.fcn = function(n, J, tau, pi00, pi01, pi10, pi11)
{
  X = rbinom(n, 1, 0.2)
  M = matrix(nrow = J, ncol = n)
  Y = matrix(nrow = J, ncol = n)
  r = sample(1:4, J, replace = T, prob = c(pi00, pi01, pi10, pi11))
  mean(r==1)
  alpha = vector()
  beta = vector()
  for(i in 1:J)
  {
    if(r[i] == 1){
      alpha[i] = 0
      beta[i] = 0
      
    }else if(r[i] ==2){
      alpha[i] = 0
      beta[i] = 0.3*tau
      
    }else if(r[i] ==3){
      alpha[i] = 0.2*tau
      beta[i] = 0
      
    }else{
      alpha[i] = 0.2*tau
      beta[i] = 0.3*tau
    }
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,] + rnorm(n)
  }
  
  
 
  
  
  ##Calculating pvalues
  p1 = vector()
  p2 = vector()
  for(i in 1:J)
  {
    p1[i] = cor.test(X, M[i,])$p.value
    p2[i] = cor.test(M[i,], Y[i,])$p.value
  }
  
  ##null estimation
  input_pvalues = cbind(p1, p2)
  nullprop = null_estimation(input_pvalues)
return(c(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2))  
}

#Parameters
###Generating Data

##Dense null
pi00 = 0.6
pi01 = 0.2
pi10 = 0.2
pi11 = 0
tau = 1 ##mediation parameter
n = 250
J = 5000


null_prop = rowMeans(replicate(20, null.fcn(n, J , tau, pi00, pi01, pi10, pi11)))

##Calculating power and FDR



fdr.fcn = function(n, J, pi00, pi01, pi10, pi11, tau)
{
  X = rbinom(n, 1, 0.2)
  M = matrix(nrow = J, ncol = n)
  Y = matrix(nrow = J, ncol = n)
  r = sample(1:4, J, replace = T, prob = c(pi00, pi01, pi10, pi11))
  mean(r==1)
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  for(i in 1:J)
  {
    if(r[i] == 1){
      alpha[i] = 0
      beta[i] = 0
      
    }else if(r[i] ==2){
      alpha[i] = 0
      beta[i] = 0.3*tau
      
    }else if(r[i] ==3){
      alpha[i] = 0.2*tau
      beta[i] = 0
      
    }else{
      alpha[i] = 0.2*tau
      beta[i] = 0.3*tau
    }
    tn[i] = alpha[i]*beta[i] ==0
    tp[i] = alpha[i]*beta[i] !=0
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,] + rnorm(n)
  }
  
  
 
  
  ##Calculating pvalues
  p1 = vector()
  p2 = vector()
  for(i in 1:J)
  {
    p1[i] = cor.test(X, M[i,])$p.value
    p2[i] = cor.test(M[i,], Y[i,])$p.value
  }
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  pval2 = DACT(p1, p2, correction = "JC")
  ##null estimation
  
  nullprop = null_estimation(input_pvalues)
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
 threshhold = max(pmax[fdr_hdmt<=0.05])
rej1 = pmax <= threshhold
fdr1 = sum(rej1*tn)/max(1,sum(rej1))
pow1 = sum(rej1*tp)/sum(tp)

rej2 = pval2 <= 0.05
fdr2 = sum(rej2*tn)/max(1,sum(rej2))
pow2 = sum(rej2*tp)/sum(tp)
return(c(fdr1,fdr2,pow1, pow2))
}

##Dense alternative
pi00 = 0.4
pi01 = 0.2
pi10 = 0.2
pi11 = 0.2
tau = 1 ##mediation parameter
n = 250
J = 5000
rowMeans(replicate(5,fdr.fcn(n, J, pi00, pi01, pi10, pi11, tau)))


##sparse alternative
pi00 = 0.88
pi01 = 0.05
pi10 = 0.05
pi11 = 0.02
tau = 1 ##mediation parameter
n = 250
J = 5000
rowMeans(replicate(5,fdr.fcn(n, J, pi00, pi01, pi10, pi11, tau)))
