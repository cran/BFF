## ---- include = FALSE---------------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE
# )

## ----setup, echo = FALSE------------------------------------------------------
library(BFF)
library(BSDA)

## -----------------------------------------------------------------------------
# generate some data
n = 100
data_one = rnorm(n = n, mean = 0.2, sd = 1)
data_two = rnorm(n = n, mean = 0.1, sd = 1)

# calculate a test statistic
z_score_one = z.test(x = data_one, sigma.x = 1)$statistic
z_score_two = z.test(x = data_one, y = data_two, sigma.x = 1, sigma.y = 1)$statistic

# calculate a BFF
z_BFF_one = BFF_z_test(z_stat = z_score_one, n = 100)
z_BFF_two = BFF_z_test(z_stat = z_score_two, one_sample = FALSE, n1 = 100, n2  = 100)


## -----------------------------------------------------------------------------
# generate some data
n = 100
data_one = rnorm(n = n, mean = -0.1)
data_two = rnorm(n = n, mean = 0.1)

# calculate a test st
t_one = t.test(x = data_one)
t_two = t.test(x = data_one, y = data_two)
t_score_one = t_one$statistic
t_score_two = t_two$statistic
t_df_one = n - 1
t_df_two = 197.9

# calculate a BFF
t_BFF_one = BFF_t_test(t_stat = t_score_one, df = t_df_one, n = 100)
t_BFF_two = BFF_t_test(t_stat = t_score_two, df = t_df_two, one_sample = FALSE, n1 = 100, n2  = 100)


## -----------------------------------------------------------------------------
# generate some data - take from chisq.test example
x <- matrix(c(12, 5, 7, 7), ncol = 2)
chi2_stat = chisq.test(x)$statistic        

# calculate a BFF
chi2_BFF_one = BFF_chi2_test(chi_stat = chi2_stat, df = 1, n = 4)

## -----------------------------------------------------------------------------
# generate some data 
n = 100
p = 3
X = matrix(rnorm(n*p), nrow = n)
beta = c(1,1,0)
y = X %*% beta + rnorm(n)
model1 = lm(y ~ X)
anova_model = anova(model1)
F_stat = anova_model$`F value`[1]

# calculate a BFF
F_BFF_one = BFF_F_test(f_stat = F_stat, df1 = anova_model$Df[1], df2 = anova_model$Df[2], n = n)

