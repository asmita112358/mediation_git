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

var_fun = function(sigma_a, sigma_b, X, M)
{
  ##...computations go here...
  m = nrow(M)
  var1 = sigma_a^2/mean(X^2)
  t1 = sigma_b^2*mean(X^2)
  var2 = vector()
  for(i in 1:m)
  {
    t2 = mean(X^2)*mean(M[i,]^2) - (mean(X*M[i,]))^2
    var2[i] = t1/t2
  }
  
  
  return(c(var1, var2))  ##Where var1= scalar and var2=vector of length m 
}
##Functions for EM algorithm.
##Computing E(gamma|everythingelse)

Q = function(alpha, beta,X, M, pi, sigma_a, sigma_b, mu, theta)
{
  m = nrow(M)
  t = matrix(nrow = m, ncol = 4)
  ##Insert code here for computing sigma1_i and sigma_2i from sigma_a and sigma_b
  variances = var_fun(sigma_a, sigma_b, X, M)
  var1 = variances[1]
  var2 = variances[-1]
  for(i in 1:m)
  {
    t[i,1] = pi[1]*dmvnorm(c(alpha[i], beta[i]), c(0,0), matrix(c(var1,0,0,var2[i]), nrow = 2))
    t[i,2] = pi[2]*dmvnorm(c(alpha[i], beta[i]), c(mu,0), matrix(c(var1 + 1,0,0,var2[i]), nrow = 2))
    t[i,3] = pi[3]*dmvnorm(c(alpha[i], beta[i]), c(0,theta), matrix(c(var1,0,0,var2[i] + 1), nrow = 2))
    t[i,4] = pi[4]*dmvnorm(c(alpha[i], beta[i]), c(mu,theta), matrix(c(var1 + 1,0,0,var2[i] + 1), nrow = 2))
    
  }
  
  return(t/rowSums(t))  ##Check this
}
