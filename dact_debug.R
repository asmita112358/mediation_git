##dact_debug

###Generating data:
tau = 1
n = 250
m = 5000
pi = c(0.49, 0.21, 0.21, 0.09) ##pi0a = 0.7, pi0b = 0.7, we get each of the elements of pi by multiplying
kap = 1
psi = 2
X = rnorm(n, 3, sd = 0.75)
generate.obj = generate(m,n,pi, tau,kap, psi, X)
M = generate.obj$M
Y = generate.obj$Y
tp = generate.obj$tp
tn = generate.obj$tn


##Computing p values
p1 = vector()
p2 = vector()
for(i in 1:m)
{
  #obj = lm(X ~ -1 + M[i,])
  #summary(obj)$coefficients[,4]
  #p1[i] = cor.test(X, M[i,])
  #obj2 = lm(Y[i,] ~ -1 + M[i,])
  #p2[i] = summary(obj2)$coefficients[,4]
  p1[i] = cor.test(X, M[i,])$p.value
  p2[i] = cor.test(M[i,], Y[i,])$p.value
}
input_pvalues = cbind(p1, p2)
pmax = apply(input_pvalues, 1, max)


##Sanity check with JS mixture
nullprop = null_estimation(input_pvalues)
fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                         nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
threshhold = max(pmax[fdr_hdmt<= size])
rej1 = pmax <= threshhold
fdr1 = sum(rej1*tn)/max(1,sum(rej1))
pow1 = sum(rej1*tp)/sum(tp)


##DACT
p_a = p1
p_b = p2
pi0a = 0.7
pi0b = 0.7
p.mat = cbind(p_a,p_b)
p3 = (apply(p.mat,1,max))^2
wg1 = pi0a*(1-pi0b)
wg2 = (1-pi0a)*pi0b
wg3 = pi0a*pi0b
wg.sum = wg1 + wg2 + wg3
wg.std = c(wg1,wg2,wg3)/wg.sum
p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
p_dact = EfronCorrect(p_dact)
rej2 = p_dact <= 0.05
fdr2 = sum(rej2*tn)/max(1,sum(rej2))
pow2 = sum(rej2*tp)/sum(tp)
