#' Sensitivity analysis for observational studies
#'
#' @inheritParams extrema.os
#' @param alpha Significance level
#' @param parallel Should parallel computing be used?
#' @param B Number of Bootstrap resamples.
#' @param warm.start Warm start the refitting of propensity score model (doesn't seem to help).
#'
#' @return A (1 - alpha) confidence interval.
#'
#' @import stats parallel
#' @export
#'

bootsens.os <- function(Z, X, Y, Lambda = 1, alpha = 0.05, estimand = "att", reg.adjust = FALSE, parallel = TRUE, B = 1000, warm.start = FALSE) {

    estimand <- match.arg(estimand, "att")

    no.cores <- ifelse(parallel, detectCores(), 1)
    n <- length(Z)

    if (warm.start) {
        start <- glm(Z ~ X, family = "binomial")$coefs
    } else {
        start <- NULL
    }

    out <- mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE);
        res <- tryCatch(extrema.os(Z[s], X[s, ], Y[s], Lambda, estimand, reg.adjust, start),
                        error = function(e) {print(e)});
        res},
        mc.cores = no.cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))

}


#' @describeIn bootsens.os Obtain extrema of IPW estimator for observational studies
#'
#' @param estimand Either "ate" (average treatment effect) or "att" (average treatment effect on the treated)
#' @inheritParams extrema.md
#'
#' @return Extrema (an interval).
#'
#' @import stats
#' @export
#'
extrema.os <- function(Z, X, Y, Lambda = 1, estimand = "att", reg.adjust = FALSE, start = NULL) {

    estimand <- match.arg(estimand, "att")
    
    if (estimand == "att") {
        if (!is.null(start)) {
            mean(Y[Z == 1]) - rev(extrema.md(Z, X, Y, Lambda, "missing", reg.adjust, - start))
        } else {
            mean(Y[Z == 1]) - rev(extrema.md(Z, X, Y, Lambda, "missing", reg.adjust))
        }
    }

}


#' @describeIn bootsens.md Obtain extrema of IPW estimator for missing data
#'
#' @inheritParams get.extrema
#' @param reg.adjust Should regression adjustment (augmented IPW) be used?
#' @param start Starting values for the propensity score model to warm start \code{glm}.
#'
#' @return Extrema (an interval).
#'
#' @import stats
#' @export
#'
extrema.md <- function(Z, X, Y, Lambda = 1, estimand = "missing", reg.adjust = FALSE, start = NULL) {
    
    estimand <- match.arg(estimand, "missing")
    
    # estimate weights
    weight_x <- ebal::ebalance(Treatment = Z, X = X)$w
    
    if (reg.adjust) {
        df <- data.frame(Y = Y, X = X)
        or.model <- lm(Y ~ ., df[Z == 1, ])
        Y.fitted <- predict(or.model, df)
    } else {
    Y.fitted <- rep(0, length(Y))
    }

    out <- get.extrema(Z, X, Y - Y.fitted, Lambda, weight_x, estimand)
    out <- out + switch(estimand, all = mean(Y.fitted), missing = mean(Y.fitted[Z == 0]))
    
}


#' Fast solution to obtain extrema of IPW estimator in missing data
#' problems
#'
#' @param Z Indicator of missingness
#' @param X A matrix of covariates
#' @param Y Outcome
#' @param Lambda Sensitivity parameter (odds ratio)
#' @param fitted.prob Fitted propensity score
#' @param estimand Target estimand, either E[Y] ("all") or E[Y|Z=0] ("missing")
#'
#' @return Extrema (an interval) of the IPW estimator.
#'
#' @export
#'
get.extrema <- function(Z, X, Y, Lambda = 1, weight_x, estimand = "missing") {
    
    estimand <- match.arg(estimand, "missing")

    # only control observations
    Y <- Y[Z == 0]
    
    # order
    weight_x <- weight_x[order(-Y)]
    Y <- Y[order(-Y)]
    
    ## maximization
    num.each.low <- Y * (Lambda^(-1) * weight_x)
    num.each.up <- Y * (Lambda * weight_x)
    num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)
    
    den.each.low <- (Lambda^(-1) * weight_x)
    den.each.up <- (Lambda * weight_x)
    den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)
    
    maximum <- max(num / den)
    
    ## minimization
    num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
    den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
    minimum <- min(num / den)
    
    c(minimum, maximum)
    
}