

source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/functions_lfdr.R")

source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/dact.R")
library(DACT)
library(HDMT)
library(locfdr)






comb.fcn = function(X, M, Y, tp, tn, size = 0.05)
{
  m = nrow(Y)
  n = ncol(Y)
  ##Calculating pvalues
  p1 = vector()
  p2 = vector()
  for(i in 1:m)
  {
    p1[i] = cor.test(X, M[i,])$p.value
    p2[i] = cor.test(M[i,], Y[i,])$p.value
  }
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  pval2 = DACT::DACT(p1, p2, correction = "JC")
  ##null estimation
  
  nullprop = null_estimation(input_pvalues)
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                           nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  threshhold = max(pmax[fdr_hdmt<= size])
  rej1 = pmax <= threshhold
  fdr1 = sum(rej1*tn)/max(1,sum(rej1))
  pow1 = sum(rej1*tp)/sum(tp)
  
  rej2 = pval2 <= size
  fdr2 = sum(rej2*tn)/max(1,sum(rej2))
  pow2 = sum(rej2*tp)/sum(tp)
  
  
  sample2 = est.coeff(X, Y, M)
  alpha = sample2[,1]
  beta = sample2[,2]
  
  ##call a function here to estimate the value of pi_start using storey's method.
  ##pi_start = pi.init(p1, p2, lambda)
  pi_start = c(nullprop$alpha00, nullprop$alpha10, nullprop$alpha01, 1 - (nullprop$alpha00 + nullprop$alpha10 + nullprop$alpha01))
  obj = maximization(alpha, beta, X, Y, M, pi_start, maxiter = 1000)
  parm = obj$par
  obj_rej = lfdr(sample2, parm)
  rej3 = obj_rej$rej
  fdr3 = sum(rej3*tn)/max(1,sum(rej3))
  pow3 = sum(rej3*tp)/sum(tp)
  
  return(c(fdr1, fdr2, fdr3, pow1, pow2, pow3))
  
}

##sparse alternative
means_mat = matrix(nrow = 3, ncol = 6)
sd_mat = matrix(nrow = 3, ncol = 6)
#tau = 0.5
#tau = 1.5
tau = 2.5
n = 100
m = 1000
pi = c(0.88, 0.05, 0.05, 0.02)
v = matrix(nrow = 20, ncol = 6)
for(l in 1:20)
  {
  X = rnorm(n, 3, sd = 0.75)
  generate.obj = generate(m,n,pi, tau, X)
  M = generate.obj$M
  Y = generate.obj$Y
  tp = generate.obj$tp
  tn = generate.obj$tn
 
 v[l,]=comb.fcn(X, M, Y, tp, tn)
 print(l)
}
means = colMeans(v, na.rm = T)
std = sqrt(apply(v,2, var, na.rm = T))


means_mat[3,] = means
sd_mat[3,] = std

write.csv(rbind(means_mat, sd_mat), file = "sparse_alt.csv")
