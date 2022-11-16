##EM function when alpha and beta are from normal distribution.

##Libraries
rm(list = ls())
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
  vec1 = rnorm(m, 0.2*tau, 1/sqrt(n))
  vec2 = rnorm(m, 0.3*tau, 1/sqrt(n))
  for(i in 1:m)
  {
    if(gamma[i] == 1){  ##h00
      alpha[i] = 0
      beta[i] = 0
      
    }else if(gamma[i] ==3){  ##h01
      alpha[i] = 0
      beta[i] = vec2[i]
      
    }else if(gamma[i] ==2){  ##h10
      alpha[i] = vec1[i]
      beta[i] = 0
      
    }else{    ##h11
      alpha[i] = vec1[i]
      beta[i] = vec2[i]
    }
    tn[i] = alpha[i]*beta[i] ==0
    tp[i] = alpha[i]*beta[i] !=0
    M[i,] = alpha[i]*X + rnorm(n)
    Y[i,] = beta[i]*M[i,] + rnorm(n)
    
  }
  return(list(M = M, Y = Y, X = X, sample = cbind(alpha, beta) ,tp = tp, tn = tn))
  
  
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
    obj2 = lm(Y[i,] ~ -1 + M[i,] )
    alpha[i] = obj1$coefficients
    beta[i] = obj2$coefficients[1]
  }
  return(cbind(alpha, beta))
}

var_fun = function(sigma_a, sigma_b, X, M)
{
  ##...computations go here...
  m = nrow(M)
  var1 = sigma_a^2/sum(X^2)
  # X = matrix(X, nrow = n)
  P = diag(n) - X%*%t(X)/sum(X^2)
  var2 = vector()
  for(i in 1:m)
  {
    #t2 = mean(X^2)*mean(M[i,]^2) - (mean(X*M[i,]))^2
    #var2[i] = t1/t2
    
    #mean(X^2/(X^2*M[i,]^2 - (X*M[i,])^2))
    var2[i] = sigma_b^2/(t(M[i,])%*%P%*%M[i,])
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
    t[i,2] = pi[2]*dmvnorm(c(alpha[i], beta[i]), c(mu,0), matrix(c(var1 + 1/n,0,0,var2[i]), nrow = 2))
    t[i,3] = pi[3]*dmvnorm(c(alpha[i], beta[i]), c(0,theta), matrix(c(var1,0,0,var2[i]+1 ), nrow = 2))
    t[i,4] = pi[4]*dmvnorm(c(alpha[i], beta[i]), c(mu,theta), matrix(c(var1 + 1/n,0,0,var2[i] + 1/n), nrow = 2))
    
  }
  
  return(t/rowSums(t))  ##Check this
}


LL.complete = function(sigma_a, sigma_b, alpha, beta, X, M, Qval, pi, mu, theta)
{
  m = nrow(M)
  
  ##Insert code here for computing sigma1_i and sigma_2i from sigma_a and sigma_b
  var1 = sigma_a^2/sum(X^2)
  # X = matrix(X, nrow = n)
  P = diag(n) - X%*%t(X)/sum(X^2)
  var2 = vector()
  for(i in 1:m)
  {
    #t2 = mean(X^2)*mean(M[i,]^2) - (mean(X*M[i,]))^2
    #var2[i] = t1/t2
    
    #mean(X^2/(X^2*M[i,]^2 - (X*M[i,])^2))
    var2[i] = sigma_b^2/(t(M[i,])%*%P%*%M[i,])
  }
  t = matrix(nrow = m, ncol = 4)
  for(i in 1:m)
  {
    t[i,1] = dmvnorm(c(alpha[i], beta[i]), c(0,0), matrix(c(var1 ,0,0,var2[i]), nrow = 2))
    t[i,2] = dmvnorm(c(alpha[i], beta[i]), c(mu,0), matrix(c(var1 + 1/n,0,0,var2[i]), nrow = 2))
    t[i,3] = dmvnorm(c(alpha[i], beta[i]), c(0,theta), matrix(c(var1,0,0,var2[i] + 1/n), nrow = 2))
    t[i,4] = dmvnorm(c(alpha[i], beta[i]), c(mu,theta), matrix(c(var1 + 1/n,0,0,var2[i] + 1/n), nrow = 2))
    
  }
  return(sum(Qval*log(t)) + sum(Qval*log(pi)))
}

LL.data = function(sigma_a, sigma_b, alpha, beta, X, M, pi, mu, theta)
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
    t[i,2] = pi[2]*dmvnorm(c(alpha[i], beta[i]), c(mu,0), matrix(c(var1 + 1/n,0,0,var2[i]), nrow = 2))
    t[i,3] = pi[3]*dmvnorm(c(alpha[i], beta[i]), c(0,theta), matrix(c(var1,0,0,var2[i] + 1/n), nrow = 2))
    t[i,4] = pi[4]*dmvnorm(c(alpha[i], beta[i]), c(mu,theta), matrix(c(var1 + 1/n,0,0,var2[i] + 1/n), nrow = 2))
    
  }
  return(sum(log(rowSums(t))))
}

maximization = function(alpha, beta, X, Y, M, pi_start, maxiter = 1000) 
{
  m = nrow(M)
  
  ##Define starting values of parameters
  
  mu_start = mean(tail(sort(alpha),m/2))
  theta_start = mean(tail(sort(beta),m/2))
  sigma_a_start = 1.5 #sqrt(var(alpha)*sum(X^2))
  sigma_b_start = 2
  Qval = Q(alpha, beta, X, M, pi_start, sigma_a_start, sigma_b_start, mu_start, theta_start)
  LL.start = LL.complete(sigma_a_start, sigma_b_start, alpha, beta, X, M, Qval,pi_start, mu_start, theta_start)
  
  ##Computing variances
  variances = var_fun(sigma_a_start, sigma_b_start, X, M)
  var1 = variances[1]
  var2 = variances[-1]
  
  ##EM update
  pi_new = colMeans(Qval)
  mu_new = sum(Qval[,2]*alpha + Qval[,4]*alpha)/sum(Qval[,2] + Qval[,4])
  v = (Qval[,3] + Qval[,4])/(var2 + 1/n)
  theta_new = sum(v*beta)/sum(v)
  sigma_a_new = optimize(LL.complete, interval = c(0.1,10) , sigma_b = sigma_b_start,alpha = alpha, beta = beta,
                         X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
  sigma_b_new = optimize(LL.complete, interval = c(0.1,10) , sigma_a = sigma_a_new,alpha = alpha, beta = beta,
                         X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
  
  ##Computing the value of log-likelihood with the EM updates
  LL.new = LL.complete(sigma_a_new, sigma_b_new, alpha, beta, X, M, pi_new, mu_new, theta_new)
  
  ##Set a counter, start a while loop
  counter = 1
  LL.vec = vector()
  LL.vec[counter] = LL.new
  while(LL.new - LL.start > 1e-1)
  {
    counter = counter + 1
    ##The new values from the last iteration becomes the starting values
    LL.start = LL.new
    LL.vec[counter] = LL.new
    print(LL.new)
    pi_start = pi_new
    mu_start = mu_new
    theta_start = theta_new
    sigma_a_start = sigma_a_new
    sigma_b_start = sigma_b_new
    Qval = Q(alpha, beta, X, M, pi_start, sigma_a_start, sigma_b_start, mu_start, theta_start)
    
    ##Computing variances
    variances = var_fun(sigma_a_start, sigma_b_start, X, M)
    var1 = variances[1]
    var2 = variances[-1]
    ##EM Update
    pi_new = colMeans(Qval)
    mu_new = sum(Qval[,2]*alpha + Qval[,4]*alpha)/sum(Qval[,2] + Qval[,4])
    v = (Qval[,3] + Qval[,4])/(var2 + 1/n)
    theta_new = sum(v*beta)/sum(v)
    sigma_a_new = optimize(LL.complete, interval = c(0.1,10) , sigma_b = 1,alpha = alpha, beta = beta,
                           X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
    sigma_b_new = optimize(LL.complete, interval = c(0.1,10) , sigma_a = sigma_a_new,alpha = alpha, beta = beta,
                           X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
    
    print(c(sigma_a_new, sigma_b_new))
    ##Increase counter and check the likelihood
    
    LL.new = LL.data(sigma_a_new, 1.7, alpha, beta, X, M, pi_new, mu_new, theta_new)
    #print(counter)
    
    ##If the while loop runs longer than maxiter then it's stopped forcibly.
    #if(counter > maxiter) {print(paste("LL=", LL.new));LL.new = LL.start}
  }
  return(list(likelihood = LL.vec, par = c(pi_new, mu_new, theta_new, sigma_a_new, sigma_b_new))) 
  
} 

##Function for computing cutoffs
lfdr = function(samp, parm, size = 0.05)
{
  J = nrow(samp)
  a = samp[,1]
  b = samp[,2]
  pi = parm[1:4]
  mu = parm[5]
  theta = parm[6]
  sigma_a = parm[7]
  sigma_b = parm[8]
  a = samp[,1]
  b = samp[,2]
  variances = var_fun(sigma_a, sigma_b, X, M)
  var1 = variances[1]
  var2 = variances[-1]
  
  lfdr = vector()
  for(i in 1:J)
  {
    t1 = pi[1]*dnorm(a[i],0,sqrt(var1))*dnorm(b[i],0,sqrt(var2[i])) + 
      pi[2]*dnorm(a[i], mu,sqrt(var1))*dnorm(b[i],0,sqrt(var2[i])) + 
      pi[3]*dnorm(a[i],0,sqrt(var1))*dnorm(b[i],theta,sqrt(var2[i]))
    t2 = pi[4]*dnorm(a[i], mu,sqrt(var1))*dnorm(b[i],theta,sqrt(var2[i]))
    lfdr[i] = t1/(t1 + t2)
    #rm(list = c(t1,t2))
  }
  st.lfdr<-sort(lfdr)
  k=1
  
  while(k<J && ((1/k)*sum(st.lfdr[1:k])) <= size){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<- lfdr<=lfdrk
  accept<- lfdr>lfdrk
  return(list(cutoff = lfdrk, rej = reject, k = k))
}

#y_plot = vector()
#x_plot = seq(0.1, 10, 0.1)
#for(i in 1:100)
#{
#  y_plot[i] = LL.data(x_plot[i],sigma_a = 1,alpha = alpha, beta = beta,
                      #   X = X, M = M, pi = pi_new, mu = mu_new, theta = theta_new )
#}

#LL.complete(1,1,alpha = alpha, beta = beta,
          #  X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new)


#optimize(LL.complete, interval = c(0.1,10) , sigma_a = sigma_a_new,alpha = alpha, beta = beta,
 #        X = X, M = M, Qval = Qval, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
#optimize(LL.data, interval = c(0.1,10) , sigma_a = sigma_a_new,alpha = alpha, beta = beta,
#         X = X, M = M, pi = pi_new, mu = mu_new, theta = theta_new, maximum = TRUE)$maximum
