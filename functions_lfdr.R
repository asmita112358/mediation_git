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
