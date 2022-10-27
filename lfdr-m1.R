##Implementation of 2d rejection region with a marginal local fdr
library(ggplot2)
##Generate data:
J = 1000
mu = 1
theta = -1
pi = c(0.7, 0.1, 0.1, 0.1)
var1 = 1
var2 = 3

obj = generate(J, pi, var1, var2, mu, theta)
samp = obj$sample
gamma = obj$gamma
tp = (gamma == 4)
tn = (gamma != 4)


##EM algorithm to estimate the parameters.
##Will write an EM algorithm to estimate the parameters later on, 
#for now, using the actual parameter values to compute the oracle rule.


lfdr = function(samp, parm)
{
  lfdra = vector()
  lfdrb = vector()
  for(i in 1:J)
  {
    t1 = (pi[1] + pi[3])*dnorm(a[i], mean = 0, sd = sqrt(var1))
    t2 = (pi[2] + pi[4])*dnorm(a[i], mean = mu, sd = sqrt(var1 + 1))
    lfdra[i] = t1/(t1 + t2)
    
    t3 = (pi[1] + pi[2])*dnorm(b[i], mean = 0, sd = sqrt(var2))
    t4 = (pi[3] + pi[4])*dnorm(b[i], mean = theta, sd = sqrt(var2 + 1))
    lfdrb[i] = t3/(t3 + t4)
  }
}

##Making a ggplot
hyp = if_else(tp, "non_null", "null")
data = data.frame(lfdra, lfdrb, hyp)
ggplot(data, aes(x = lfdra, y = lfdrb, color = hyp)) + geom_point()
