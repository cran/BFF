## ----include = FALSE----------------------------------------------------------
# knitr::opts_chunk$set(
#   collapse = TRUE
# )

## ----setup, echo = FALSE------------------------------------------------------
library(BFF)
library(BSDA)

## -----------------------------------------------------------------------------
# Define the density function
density_function <- function(lambda, tau2, r) {
  coef <- (lambda^2)^r / ((2 * tau2)^(r + 0.5) * gamma(r + 0.5))
  exp_term <- exp(-lambda^2 / (2 * tau2))
  return(coef * exp_term)
}

# Define omega and n
omega <- 0.5
n <- 50

# Define the range for lambda
lambda_range <- seq(-10, 10, length.out = 1000)

# Define the values of r to plot
r_values <- c(1, 2, 3, 5, 10)

# Define colors for the different r values
colors <- rainbow(length(r_values))

# Plot the densities using base R plotting functions
plot(NULL, xlim = c(-10, 10), xlab = expression(lambda), ylab = "Density", main = "Density for varying r", ylim = c(0, 0.4))

# Add lines for each value of r
for (i in seq_along(r_values)) {
  r <- r_values[i]
  tau2 <- (n * omega^2) / (2 * r)
  densities <- sapply(lambda_range, density_function, tau2 = tau2, r = r)
  lines(lambda_range, densities, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = paste("r =", r_values), col = colors, lwd = 2)


## -----------------------------------------------------------------------------
tBFF = t_test_BFF(t_stat = 2.5, n = 50, one_sample = TRUE)
tBFF # plot the results
plot(tBFF)

# now an example with different types of settings
tBFF_twosample = t_test_BFF(t_stat = 2.5, n1 = 50, n2 = 14, one_sample = FALSE)
tBFF_twosample # plot the results
plot(tBFF_twosample)


# repeat the above with r = 5
tBFF_r5 = t_test_BFF(t_stat = 2.5, n = 50, one_sample = TRUE, r = 5)
tBFF_r5 # plot the results
plot(tBFF_r5)

# now an example with different types of settings
tBFF_twosample_r5 = t_test_BFF(t_stat = 2.5, n1 = 50, n2 = 14, one_sample = FALSE, r = 5)
tBFF_twosample_r5 # plot the results
plot(tBFF_twosample_r5)

