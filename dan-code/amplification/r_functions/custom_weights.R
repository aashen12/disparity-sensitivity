################################################################
## Code to do bootstrap sensitivity with custom weight functions
################################################################

#' Fast solution to obtain extrema of IPW estimator in missing data
#' problems
#'
#' @param A Indicator of missingness
#' @param X A matrix of covariates
#' @param Y Outcome
#' @param gamma Sensitivity parameter (log odds ratio)
#' @param fitted.prob Fitted propensity score
#' @param estimand Target estimand, either E[Y] ("all") or E[Y|A=0] ("missing")
#'
#' @return Extrema (an interval) of the IPW estimator.
#'
#' @export
#'
get_extrema <- function(A, X, Y, weight_x, gamma = 0) {

    
    eg <- weight_x
    # CHANGE A == 0
    Y <- Y[A == 0]
    # eg already only just control observations
    # eg <- eg[A == 0]
    eg <- eg[order(-Y)]
    Y <- Y[order(-Y)]

    ## maximization
    # CHANGE: remove c
    num.each.low <- Y * (exp(-gamma) * eg)
    num.each.up <- Y * (exp(gamma) * eg)
    num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

    den.each.low <- (exp(-gamma) * eg)
    den.each.up <- (exp(gamma) * eg)
    den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

    maximum <- max(num / den)
    ## print(den[which.max(num/den)] / n)

    ## minimization
    num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
    den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
    minimum <- min(num / den)
    ## print(den[which.min(num/den)] / n)

    c(minimum, maximum)

}


extrema_md_weights <- function(A, X, Y, weightfunc,
                               gamma = 0) {


  weight_x <- weightfunc(A, X)

  Y.fitted <- rep(0, length(Y))

# CHANGE: weight_x
  out <- get_extrema(A, X, Y - Y.fitted, weight_x, gamma)
  out <- out + mean(Y.fitted[A == 0])
  return(out)
}


#' @describeIn bootsens.os Obtain extrema of IPW estimator for observational studies
#'
#' @param estimand Either "ate" (average treatment effect) or "att" (average treatment effect on the treated)
#' @inheritParams extrema.md
#'
#' @return Extrema (an interval).
#'
#' @import stats
extrema_os <- function(A, X, Y, weightfunc, gamma = 0) {


  mean(Y[A == 1]) - rev(extrema_md_weights(A, X, Y, weightfunc, gamma))

}


#' Sensitivity analysis for observational studies
#'
#' @inheritParams extrema.os
#' @param alpha Significance level
#' @param B Number of Bootstrap resamples.
#' @param warm.start Warm start the refitting of propensity score model (doesn't seem to help).
#'
#' @return A (1 - alpha) confidence interval.
bootsens_os <- function(A, X, Y, weightfunc, gamma = 0, alpha = 0.05, n_cores = 1, B = 1000) {

    n <- length(A)

    out <- parallel::mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE);
        res <- tryCatch(extrema_os(A[s], X[s, ], Y[s], weightfunc, gamma),
                        error = function(e) {print(e)});
        res},
        mc.cores = n_cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))

}
