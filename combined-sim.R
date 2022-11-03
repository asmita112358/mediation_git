##Function for generating from the original model

generate = function(m, n, pi, tau, X)
{
  pi00 = pi[1]
  pi10 = pi[2]
  pi01 = pi[3]
  pi11 = pi[4]
  M = matrix(nrow = J, ncol = n)
  Y = matrix(nrow = J, ncol = n)
  gamma = sample(1:4, J, replace = T, prob = c(pi00, pi01, pi10, pi11))
  
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
  
  
}