
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/dact.R")
source("~/Documents/OneDrive - Texas A&M University/Documents/fair-ML/RCodes/mediation_git/EM2v2.R")



n = 100
m = 1000
pi = c(0.7, 0.1, 0.1, 0.1)
tau = 2
X = rnorm(n, 3, sd = 0.75)
generate.obj = generate(m,n,pi, tau, X)
attach(generate.obj)


comb.fcn = function(X, M, Y, tp, tn, size)
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
  pval2 = DACT(p1, p2, correction = "JC")
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
  
  
  
}
