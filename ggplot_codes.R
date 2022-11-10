library(ggplot2)
library(dplyr)
library(tidyr)
##Setup 1: Known, fixed alpha, beta, under dense and sparse slternative

data = read.csv("sparse_alt.csv")
method = c("JS-mixture","DACT","MLFDR")
tau = c("low", "medium", "high")
fdr = data[ 1:3, 2:4]
fdr_se = data[4:6, 2:4]
pow = data[1:3, 5:7]
pow_se = data[4:6, 5:7]


colnames(fdr) = method
rownames(fdr) = tau
fdr_matrix_long = gather(fdr, key = "method", value = "fdr")
fdr_matrix_long$tau = factor(rep(tau,3), levels = c("low", "medium", "high"))
se = as.vector(as.matrix(fdr_se))
fdr_matrix_long$se = se


ggplot(fdr_matrix_long, aes(x = tau, y = fdr, fill = method))+
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0.05)+
  geom_errorbar(aes(ymin = fdr - 1.96*se/sqrt(20), ymax = fdr + 1.96*se/sqrt(20)), width = 0.2,position = position_dodge(width = 0.9)) +
  labs(y = "FDR", x = "Degree of Mediation")
  

fdr = pow
fdr_se = pow_se
colnames(fdr) = method
rownames(fdr) = tau
fdr_matrix_long = gather(fdr, key = "method", value = "fdr")
fdr_matrix_long$tau = factor(rep(tau,3), levels = c("low", "medium", "high"))
se = as.vector(as.matrix(fdr_se))
fdr_matrix_long$se = se


ggplot(fdr_matrix_long, aes(x = tau, y = fdr, fill = method))+
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = fdr - 1.96*se/sqrt(20), ymax = fdr + 1.96*se/sqrt(20)), width = 0.2,position = position_dodge(width = 0.9)) +
  labs(y = "Power", x = "Degree of Mediation")
  