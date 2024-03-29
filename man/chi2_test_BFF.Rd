% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FINAL_FUNCTIONS_chi2_test.R
\name{chi2_test_BFF}
\alias{chi2_test_BFF}
\title{chi2_test_BFF}
\usage{
chi2_test_BFF(
  chi2_stat,
  n = NULL,
  df = NULL,
  pearsons = TRUE,
  savename = NULL,
  maximize = FALSE,
  r = 1,
  tau2 = NULL,
  save = TRUE,
  xlab = NULL,
  ylab = NULL,
  main = NULL
)
}
\arguments{
\item{chi2_stat}{chi^2 statistic}

\item{n}{sample size}

\item{df}{degrees of freedom}

\item{pearsons}{Is this a test of Pearson’s chi^2 test for goodness-of-fit? Default is TRUE. FALSE assumes a likelihood ratio test}

\item{savename}{optional, filename for saving the pdf of the final plot}

\item{maximize}{Should the value of r be maximized? Default is FALSE. Only set to TRUE if analyzing multiple studies}

\item{r}{r value}

\item{tau2}{tau2 values (can be a single entry or a vector of values)}

\item{save}{should a copy of the plot be saved?}

\item{xlab}{optional, x label for plot}

\item{ylab}{optional, y label for plot}

\item{main}{optional, main label for plot}
}
\value{
Returns Bayes factor function results
 \tabular{ll}{
   \code{BFF} \tab The log of the Bayes Factor Function values \cr
   \tab \cr
   \code{effect_size} \tab Effect sizes tested (seq(0, 1, by = 0.01)) \cr
   \tab \cr
   \code{BFF_max_RMSE} \tab Maximum BFF value \cr
   \tab \cr
   \code{max_RMSE} \tab Effect size that maximizes BFF\cr
   \tab \cr
   \code{tau2} \tab tau^2 values tested\cr
}
}
\description{
chi2_test_BFF constructs BFFs based on the chi-squared test. BFFs depend on hyperparameters r and tau^2 which determine the shape and scale of the prior distributions which define the alternative hypotheses.
By setting r > 1, we use higher-order moments for replicated studies. Fractional moments are set with r > 1 and r not an integer.
All results are on the log scale.
Plot saved to working directory unless a full path is specified in the 'savename' variable of the function.
}
\examples{
chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE)
chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = 0.5)
chi2BFF = chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, save = FALSE, tau2 = c(0.5, 0.8))
chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, pearsons = FALSE, save = FALSE)
chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2, save = FALSE)
chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2, pearsons = FALSE, save = FALSE)
chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2.5, save = FALSE)
chi2_test_BFF(chi2_stat = 2.5, n = 50, df = 49, r = 2.5, pearsons = FALSE, save = FALSE)
chi2_test_BFF(chi2_stat=2.5, n = 50, df = 49, maximize = TRUE)
chi2_test_BFF(chi2_stat=2.5, n = 50,  df = 49, maximize = TRUE, tau2 = 0.5)
chi2_test_BFF(chi2_stat=2.5, n = 50,  df = 49, maximize = TRUE, tau2 = c(0.5, 0.8))
chi2BFF$BFF_max_RMSE  # maximum BFF value
chi2BFF$max_RMSE      # effect size which maximizes the BFF
}
