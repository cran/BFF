#' @title Summarize BFF object
#'
#' @param object a BFF object
#' @param ... additional arguments (unused)
#'
#' @return prints summary of a BFF object.
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#' @export
summary.BFF <- function(object, ...){
  print.BFF(object, ...)
}

#' @title Summarize BFF object
#'
#' @param x a BFF object
#' @param ... additional arguments (unused)
#'
#' @return prints summary of a BFF object.
#'
#' @seealso [z_test_BFF()], [t_test_BFF()], [chi2_test_BFF()], [f_test_BFF()]
#' @export
print.BFF <- function(x, ...) {
  cat(paste0("\n\t", .test_type_name(x$test_type, x$input$one_sample)))
  cat("\n\n")

  # # print in favor of alternative
  # cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$log_bf_h1))
  # if(x$generic_test){
  #   cat(gettextf("%1$s tau2 = %2$.2f\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$omega_h1))
  # }else{
  #   cat(gettextf("%1$somega = %2$.2f (%3$s)\n", if(!x$omega_set) "maximized (in favor of alternative) " else "", x$omega_h1, .test_effect_size_name(x$test_type)))
  # }
  #
  # # print in favor of null, only if no omega is set
  # if (!x$omega_set) {
  #   cat(gettextf("%1$slog Bayes factor = %2$.2f\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$log_bf_h0))
  #   if(x$generic_test){
  #     cat(gettextf("%1$s tau2 = %2$.2f\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$omega_h0))
  #   }else{
  #     cat(gettextf("%1$somega = %2$.2f (%3$s)\n", if(!x$omega_set) "minimized (in favor of null for medium/large effect sizes) " else "", x$omega_h0, .test_effect_size_name(x$test_type)))
  #   }
  # }


  # if(!is.null(x$input$alternative)) cat(paste0("alternative = ", x$input$alternative.original))
  ## checking if there is switching
  # diff(sign(x$BFF$log_bf)) # is -2 when transitioning from positive to negative.


  # printing when omega is not set
  if (!x$omega_set) {
    has_switched = any(diff(sign(x$BFF$log_bf)) == -2)
    indices = which(diff(sign(x$BFF$log_bf)) == -2) + 1
    if (all(x$BFF$log_bf <= 0)) {
      cat("\n")
      cat("The BF provides evidence for the null hypothesis across all standardized effect sizes")
    } else {
      cat(paste0("The log BF is maximized in favor of the alternative hypothesis at the standardized effect size of ", x$omega_h1, " with value ", round(x$log_bf_h1,2), ". Standardized effect size should be chosen for individual hypotheses based on scientific intent and plausibilty."))
      cat("\n")
      if(has_switched) {
        cat("\n")
        cat(paste("The BF switches from providing evidence for the alternative hypothesis to evidence for the null hypothesis at the standardized effect size of", x$BFF$omega[indices]))
      }
    }

  } else {
    cat(paste("Log BF =", round(x$log_bf_h1, 2), "at standardized effect size", x$omega_h1))
  }

  cat("\n\n")
  # if(!is.null(x$input$alternative)) cat(paste0("alternative = ", x$input$alternative.original))

}

.test_type_name <- function(test_type, one_sample) {
  starting_strng = gettextf("Bayes Factor Test With Non-Local Priors For A %1$s%2$s",
                            if(!is.null(one_sample)) {if(one_sample) "One-Sample " else "Two-Sample "} else "",
                            switch(test_type,
                                   "t_test"           = "t-Test",
                                   "z_test"           = "z-Test",
                                   "chi2_test"        = "chi2-Test",
                                   "f_test"           = "F-Test",
                                   "regression_test"  = "Regression Test",
                                   "correlation_test" = "Correlation Test"))

}
.test_effect_size_name <- function(test_type){
  switch(test_type,
         "t_test"           = "Cohen's d",
         "z_test"           = "Cohen's d",
         "f_test"           = "Cohen's f",
         "chi2_test"        = "Cohen's w",
         "regression_test"  = "Cohen's d",
         "correlation_test" = "correlation coefficient")
}

.test_effect_size_name_plot <- function(test_type){
  switch(test_type,
         "t_test"           = "|Cohen's d|",
         "z_test"           = "|Cohen's d|",
         "f_test"           = "|Cohen's f|",
         "chi2_test"        = "|Cohen's w|",
         "regression_test"  = "|Cohen's d|",
         "correlation_test" = "correlation coefficient")
}

