# mediation_git

EM2v2.R contains the functions required to 
a) generate from a mixture normal
b) run and EM algorithm and estimate the parameters associated with the mixture normal.

lfdr-m1 contains the two dimensional rejection region based on the marginal lfdrs
lfdr-m2 contains the one dimensional rejection region based on the joint lfdr.

lfdra-vs-b contains a plot showing the marginal lfdrs for a and b for pi = c(0.7, 0.1, 0.1, 0.1), mu = 1, theta = -1, var1 = 1, var2 = 3. 
Clearly, for this setup, the margina lfdrs are unable to provide a distinction between the nulls and the non nulls.
