# mediation_git

1. EM2v2.R contains the functions required to 
 a) generate from a mixture normal
 b) run and EM algorithm and estimate the parameters associated with the mixture normal.

2. lfdr-m1 contains the two dimensional rejection region based on the marginal lfdrs.
3. lfdr-m2 contains the one dimensional rejection region based on the joint lfdr.


4. lfdra-vs-b contains a plot showing the marginal lfdrs for a and b for pi = c(0.7, 0.1, 0.1, 0.1), mu = 1, theta = -1, var1 = 1, var2 = 3. 
Clearly, for this setup, the margina lfdrs are unable to provide a distinction between the nulls and the non nulls.

5. lfdra-vs-b-setup2 contains the plot showing the marginal lfdrs for a and b for pi = (0.7,0.1, 0.1, 0.1) and mu = 3, theta = -3, var1 = 1, var2 = 2. This probably indicates that lfdr will be able to distinguish between the null and the alternative more effectively if the means of the non-null cases are far away from zero. 

6. The estimate of the mFDR and mFNR in the 2d method is incorrect. (For details, check overleaf.)Edit: this has been corrected. There was an issue with 1/m.

7. (11/7) Currently the EM function still has issues with estimating sigma_a and sigma_b. Update: This issue has been fixed.
