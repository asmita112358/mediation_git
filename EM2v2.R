##EM2_fresh

rm(list = ls())
library(MASS)
library(dplyr)
library(emdbook)
library(philentropy)
library(rootSolve)
library(dplyr)
require(EnvStats)
library(stats)
###Function to generate from a mixture normal distribution

generate = function(J, pi, var1, var2, mu, theta )
{
  tau = rnorm(J, mean = mu)
  psi = rnorm(J, mean = theta)
  gamma = sample(1:4,J, replace = TRUE, prob = pi)
  sample = matrix(nrow = J, ncol = 2)
  for(i in 1:J)
  {
    if(gamma[i] == 1){
      sample[i,] = mvrnorm(1, c(0,0),  matrix(c(var1,0,0,var2), nrow = 2))
    }else if(gamma[i] == 2){
      sample[i,] = mvrnorm(1, c(mu,0), matrix(c(var1 + 1,0,0,var2), nrow = 2))
    }else if(gamma[i] == 3){
      sample[i,] = mvrnorm(1, c(0,theta), matrix(c(var1,0,0,var2 + 1), nrow = 2))
    }else{
      sample[i,] = mvrnorm(1, c(mu,theta), matrix(c(var1 + 1,0,0,var2 + 1), nrow = 2))
    }
  }
  return(list(gamma = gamma, sample = sample))
}


##Computing E(gamma|everythingelse)

Q = function(sample, pi, var1, var2, mu, theta)
{
  t = matrix(nrow = J, ncol = 4)
  t[,1] = pi[1]*dmvnorm(sample, c(0,0), matrix(c(var1,0,0,var2), nrow = 2))
  t[,2] = pi[2]*dmvnorm(sample, c(mu,0), matrix(c(var1 + 1,0,0,var2), nrow = 2))
  t[,3] = pi[3]*dmvnorm(sample, c(0,theta), matrix(c(var1,0,0,var2 + 1), nrow = 2))
  t[,4] = pi[4]*dmvnorm(sample, c(mu,theta), matrix(c(var1 + 1,0,0,var2 + 1), nrow = 2))
  return(t/rowSums(t))
}




LL = function(var1, var2, Qval,sample, pi, mu, theta)
{
  t = matrix(nrow = J, ncol = 4)
  t[,1] = dmvnorm(sample, c(0,0), matrix(c(var1,0,0,var2), nrow = 2))
  t[,2] = dmvnorm(sample, c(mu,0), matrix(c(var1 + 1,0,0,var2), nrow = 2))
  t[,3] = dmvnorm(sample, c(0,theta), matrix(c(var1,0,0,var2 + 1), nrow = 2))
  t[,4] = dmvnorm(sample, c(mu,theta), matrix(c(var1 + 1,0,0,var2 + 1), nrow = 2))
  return(sum(Qval*log(t)) + sum(Qval*log(pi)))
}

LL.data = function(var1, var2, sample, pi, mu, theta)
{
  t = matrix(nrow = J, ncol = 4)
  t[,1] = dmvnorm(sample, c(0,0), matrix(c(var1,0,0,var2), nrow = 2))
  t[,2] = dmvnorm(sample, c(mu,0), matrix(c(var1 + 1,0,0,var2), nrow = 2))
  t[,3] = dmvnorm(sample, c(0,theta), matrix(c(var1,0,0,var2 + 1), nrow = 2))
  t[,4] = dmvnorm(sample, c(mu,theta), matrix(c(var1 + 1,0,0,var2 + 1), nrow = 2))
  return(sum(log(rowSums(pi*t))))
}

maximization = function(sample, pi_start, maxiter = 1000) 
{
  J = nrow(sample)
  alpha = sample[,1]
  beta = sample[,2]
  
  ##Define starting values of parameters
  
  mu_start = mean(alpha)
  theta_start = mean(beta)
  var1_start = var(alpha)
  var2_start = var(beta)
  Qval = Q(sample, pi_start, var1_start, var2_start, mu_start, theta_start)
  LL.start = LL.data(var1_start, var2_start, sample, pi_start, mu_start, theta_start)
  
  ##EM update
  pi_new = colMeans(Qval)
  mu_new = sum(Qval[,2]*alpha + Qval[,4]*alpha)/sum(Qval[,2] + Qval[,4])
  theta_new = sum(Qval[,3]*beta + Qval[,4]*beta)/sum(Qval[,3] + Qval[,4])
  var1_new = optimize(LL, interval = c(0,10) , var2 = 1,Qval = Qval,sample = sample, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
  var2_new = optimize(LL, interval = c(0,10) , var1 = 1,Qval = Qval,sample = sample, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
  
  ##Computing the value of log-likelihood with the EM updates
  LL.new = LL.data(var1_new, var2_new, sample, pi_new, mu_new, theta_new)
  
  ##Set a counter, start a while loop
  counter = 1
  LL.vec = vector()
  LL.vec[counter] = LL.new
  while(LL.new - LL.start > 1e-6)
  {
    counter = counter + 1
    ##The new values from the last iteration becomes the starting values
    LL.start = LL.new
    LL.vec[counter] = LL.new
    pi_start = pi_new
    mu_start = mu_new
    theta_start = theta_new
    var1_start = var1_new
    var2_start = var2_new
    Qval = Q(sample, pi_start, var1_start, var2_start, mu_start, theta_start)
    
    
    ##EM Update
    pi_new = colMeans(Qval)
    mu_new = sum(Qval[,2]*alpha + Qval[,4]*alpha)/sum(Qval[,2] + Qval[,4])
    theta_new = sum(Qval[,3]*beta + Qval[,4]*beta)/sum(Qval[,3] + Qval[,4])
    var1_new = optimize(LL, interval = c(0,10) , var2 = 1,Qval = Qval,sample = sample, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
    var2_new = optimize(LL, interval = c(0,10) , var1 = 1,Qval = Qval,sample = sample, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
    
    
    ##Increase counter and check the likelihood
   # counter = counter + 1
    LL.new = LL.data(var1_new, var2_new, sample, pi_new, mu_new, theta_new)
    #print(counter)
    
    ##If the while loop runs longer than maxiter then it's stopped forcibly.
    if(counter > maxiter) {print(paste("LL=", LL.new));LL.new = LL.start}
  }
  return(list(likelihood = LL.vec, par = c(pi_new, mu_new, theta_new, var1_new, var2_new))) 
  
}  




