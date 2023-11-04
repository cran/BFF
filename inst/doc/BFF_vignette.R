## ---- include = FALSE---------------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE
# )

## ----setup, echo = FALSE------------------------------------------------------
library(BFF)
library(BSDA)

## ----z-statistics-------------------------------------------------------------
# generating some data
n = 100
data_one = rnorm(n = n, mean = 0.2, sd = 1)
data_two = rnorm(n = n, mean = 0.1, sd = 1)

# calculating test statistics using z.test
# one-sample z-test
z_score_one = z.test(x = data_one, sigma.x = 1)$statistic
# two-sample z-test
z_score_two = z.test(x = data_one, y = data_two, sigma.x = 1, sigma.y = 1)$statistic


## ----calculating BFF for z test-----------------------------------------------
# default r and tau2
z_BFF_one = z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE) #one sample z-test
z_BFF_two = z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE) #two sample z-test
  
# default r and user specified tau2
# single tau2
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, tau2 = 0.5) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, tau2 = 0.5) #two sample z-test
# vector of tau2 values
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, tau2 = c(0.5, 0.8)) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, tau2 = c(0.5, 0.8)) #two sample z-test
  
# user specified r and default tau2
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, r = 2) #one sample z-test, integer r >1  (higher order moments) 
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, r = 2) #two sample z-test, integer r >1  (higher order moments) 
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, r = 2.5) #one sample z-test, continuous r (fractional moments)
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, r = 2.5) #two sample z-test, continuous r (fractional moments)

## ----t-statistics-------------------------------------------------------------
# generating some data
n = 100
data_one = rnorm(n = n, mean = -0.1)
data_two = rnorm(n = n, mean = 0.1)

# calculating test statistics using t.test
t_one = t.test(x = data_one)
t_two = t.test(x = data_one, y = data_two)
t_score_one = t_one$statistic
t_score_two = t_two$statistic
t_df_one = n - 1
t_df_two = 197.9


## ----calculating BFF for t test-----------------------------------------------
# default r and tau2
t_BFF_one = t_test_BFF(t_stat = t_score_one, df = t_df_one, n = 100, save = FALSE) #one sample t-test
t_BFF_two = t_test_BFF(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100, save = FALSE) #two sample t-test
  
# default r and user specified tau2
# single tau2
t_test_BFF(t_stat = t_score_one, df = t_df_one, n = 100, save = FALSE, tau2 = 0.5) #one sample t-test
t_test_BFF(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100, save = FALSE, tau2 = 0.5) #two sample t-test 
# vector of tau2 values
t_test_BFF(t_stat = t_score_one, df = t_df_one, n = 100, save = FALSE, tau2 = c(0.5, 0.8)) #one sample t-test
t_test_BFF(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100, save = FALSE, tau2 = c(0.5, 0.8)) #two sample t-test 
  
# user specified r and default tau2
t_test_BFF(t_stat = t_score_one, df = t_df_one, n = 100, save = FALSE, r = 2) #one sample t-test, integer r >1  (higher order moments) 
t_test_BFF(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100, save = FALSE, r = 2) #two sample t-test, integer r >1  (higher order moments) 
t_test_BFF(t_stat = t_score_one, df = t_df_one, n = 100, save = FALSE, r = 2.5) #one sample t-test, continuous r (fractional moments)
t_test_BFF(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100, save = FALSE, r = 2.5) #two sample t-test, continuous r (fractional moments)

## ----chi2-test----------------------------------------------------------------
# generate some data
x <- matrix(c(12, 5, 7, 7), ncol = 2)

# calculating chi2 test statistic from chisq.test
chi2_stat = chisq.test(x)$statistic        

## ----calculating BFF for chi2 test--------------------------------------------
# default r and tau2
chi2_BFF_pear = chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE) #Pearson's chi2 test
chi2_BFF_lrt = chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, pearsons = FALSE) #Likelihood ratio chi2 test
  
# default r and user specified tau2
# single tau2
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, tau2 = 0.5) #Pearson's chi2 test
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, pearsons = FALSE, tau2 = 0.5) #Likelihood ratio chi2 test
# vector of tau2 values
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, tau2 = c(0.5, 0.8)) #Pearson's chi2 test
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, pearsons = FALSE, tau2 = c(0.5, 0.8)) #Likelihood ratio chi2 test
  
# user specified r and default tau2
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, r = 2) #Pearson's chi2 test, integer r >1  (higher order moments) 
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, pearsons = FALSE, r = 2) #Likelihood ratio chi2 test, integer r >1  (higher order moments) 
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, r = 2.5) #Pearson's chi2 test, continuous r (fractional moments) 
chi2_test_BFF(chi2_stat = chi2_stat, df = 1, n = 4, save = FALSE, pearsons = FALSE, r = 2.5) #Likelihood ratio chi2 test, continuous r (fractional moments) 

## ----f statistics-------------------------------------------------------------
# generate some data
n = 100
p = 3
X = matrix(rnorm(n*p), nrow = n)
beta = c(1,1,0)
y = X %*% beta + rnorm(n)
model1 = lm(y ~ X)
anova_model = anova(model1)
F_stat = anova_model$`F value`[1]

## ----calculating BFF for f test-----------------------------------------------
# default r and tau2
F_BFF_one = f_test_BFF(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n, save = FALSE)

# default r and user specified tau2
# single tau2
f_test_BFF(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n, tau2 = 0.5, save = FALSE)
# vector of tau2 values
f_test_BFF(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n, tau2 = c(0.5, 0.8), save = FALSE)
  
# user specified r and default tau2
f_test_BFF(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n, r = 2, save = FALSE) #integer r >1  (higher order moments) 
f_test_BFF(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n, r = 2.5, save = FALSE) #continuous r (fractional moments)

## ----maximing r for z test----------------------------------------------------
# default tau2
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, maximize = TRUE) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, maximize = TRUE) #two sample z-test

# user specified tau2
#single tau2
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, tau2 = 0.5, maximize = TRUE) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, tau2 = 0.5, maximize = TRUE) #two sample z- test
# vector of tau2 values
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, tau2 = c(0.5, 0.8), maximize = TRUE) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = FALSE, tau2 = c(0.5, 0.8), maximize = TRUE) #two sample z-test

## ----plotting for z test------------------------------------------------------
# saving the plot as a pdf with default name (BFF_plot.pdf). Stored in working directory.
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE) #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = TRUE) #two sample z-test

# saving the plot as a pdf with user specified name. 
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, savename = "z-BFF-one.pdf") #one sample z-test
z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = TRUE, savename = "z-BFF-two.pdf") #two sample z-test
 
# customizing x-axis labels, y-axis labels and main title
z_test_BFF(z_stat = z_score_one, n = 100, save = FALSE, xlab = "RMSE", ylab = "Logarithm of Bayes Factor", main = "BFF curves") #one sample z-test
z_BFF_two = z_test_BFF(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2 = 100, save = TRUE, xlab = expression(tilde(omega)), ylab = expression(log(BF[10])), main = "BFF curves") #two sample z-test

