##Libraries
library(MASS)
library(dplyr)
library(emdbook)
library(philentropy)
library(rootSolve)
library(dplyr)
require(EnvStats)
library(stats)

##Function for generating from the original model
generate = function(m, n, pi, tau, X)
{
  if(length(X) != n)stop(print("Length of X and n is not same"))
  M = matrix(nrow = m, ncol = n)
  Y = matrix(nrow = m, ncol = n)
  gamma = sample(1:4, m, replace = T, prob = pi)
  
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  for(i in 1:m)
  {
    if(gamma[i] == 1){  ##h00
      alpha[i] = 0
      beta[i] = 0
      
    }else if(gamma[i] ==3){  ##h01
      alpha[i] = 0
      beta[i] = 0.3*tau
      
    }else if(gamma[i] ==2){  ##h10
      alpha[i] = 0.2*tau
      beta[i] = 0
      
    }else{    ##h11
      alpha[i] = 0.2*tau
      beta[i] = 0.3*tau
    }
    tn[i] = alpha[i]*beta[i] ==0
    tp[i] = alpha[i]*beta[i] !=0
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,] + rnorm(n)
    
  }
  return(list(M = M, Y = Y, X = X, tp = tp, tn = tn))
  
  
}
##function for estimating alpha and beta by fitting regression equations

est.coeff = function(X, Y, M) ##Takes in the original sample, returns alpha and beta
{
  m = nrow(Y)
  n = ncol(Y)
  alpha = vector()
  beta = vector()
  for(i in 1:m)
  {
    obj1 = lm(M[i,] ~ -1 + X)
    obj2 = lm(Y[i,] ~ -1 + M[i,])
    alpha[i] = obj1$coefficients
    beta[i] = obj2$coefficients
  }
  return(cbind(alpha, beta))
}


##Functions for EM algorithm.
##Computing E(gamma|everythingelse)

Q = function(sample, pi, sigma_a, sigma_b, mu, theta)
{
  t = matrix(nrow = J, ncol = 4)
  ##Insert code here for computing sigma1_i and sigma_2i from sigma_a and sigma_b
  
  
  t[,1] = pi[1]*dmvnorm(sample, c(0,0), matrix(c(var1,0,0,var2), nrow = 2))
  t[,2] = pi[2]*dmvnorm(sample, c(mu,0), matrix(c(var1 + 1,0,0,var2), nrow = 2))
  t[,3] = pi[3]*dmvnorm(sample, c(0,theta), matrix(c(var1,0,0,var2 + 1), nrow = 2))
  t[,4] = pi[4]*dmvnorm(sample, c(mu,theta), matrix(c(var1 + 1,0,0,var2 + 1), nrow = 2))
  
  return(t/rowSums(t))
}
