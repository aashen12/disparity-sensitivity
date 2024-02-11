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
get.extrema <- function(A, X, Y, gamma = 0, weight_x, estimand = c("all", "missing")) {
#CHANGE weight_x
    
    #estimand <- match.arg(estimand, c("all", "missing"))
    #c <- as.numeric(estimand == "all")

    #fitted.logit <- qlogis(fitted.prob)

    #eg <- exp(-fitted.logit)
    
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

bootsens.os <- function(A, X, Y, gamma = 0, alpha = 0.05, estimand = c("ate", "att"), reg.adjust = FALSE, parallel = TRUE, B = 1000, warm.start = FALSE) {

    estimand <- match.arg(estimand, c("ate", "att"))

    no.cores <- ifelse(parallel, detectCores(), 1)
    n <- length(A)

    if (warm.start) {
        start <- glm(A ~ X, family = "binomial")$coefs
    } else {
        start <- NULL
    }

    out <- mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE);
        res <- tryCatch(extrema.os(A[s], X[s, ], Y[s], gamma, estimand, reg.adjust, start),
                        error = function(e) {print(e)});
        res},
        mc.cores = no.cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))

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
extrema.md <- function(A, X, Y, gamma = 0, estimand = c("all", "missing"), reg.adjust = FALSE, start = NULL) {

    estimand <- match.arg(estimand, c("all", "missing"))
    # CHANGE: estimate weights instead of pscore
    # ps.model <- glm(A ~ X, family = "binomial", start = start)
    # fitted.prob <- predict(ps.model, type = "response")
    
    #weight_x <- ebalance(Treatment = A, X = X)$w
    #weight_x <- balancer::multilevel_qp(X = X, trt = A, Z = rep(1, nrow(X)), exact_global = F)$weights[A==0]
    data_frame = as.data.frame(cbind(A, X, Y))
    # Define treatment indicator and
    z_ind = "A"
    # moment covariates
    bal = list()
    bal$bal_cov = colnames(X)
    # Set tolerances
    bal$bal_tol = 0.02
    bal$bal_std = "group"
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    bal$bal_alg = FALSE
    sbwatt_object = sbw(dat = data_frame, ind = z_ind, out = "Y", bal = bal)

    weight_x <- sbwatt_object$dat_weights[sbwatt_object$dat_weights$A == 0,]$sbw_weights
    # if (reg.adjust) {
    #     df <- data.frame(Y = Y, X = X)
    #     or.model <- lm(Y ~ ., df[A == 1, ])
    #     Y.fitted <- predict(or.model, df)
    # } else {
         Y.fitted <- rep(0, length(Y))
    # }
# CHANGE: weight_x
    out <- get.extrema(A, X, Y - Y.fitted, gamma, weight_x, estimand)
    out <- out + switch(estimand, all = mean(Y.fitted), missing = mean(Y.fitted[A == 0]))

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
extrema.os <- function(A, X, Y, gamma = 0, estimand = c("ate", "att"), reg.adjust = FALSE, start = NULL) {

    estimand <- match.arg(estimand, c("ate", "att"))
# CHANGE: 1 - A to A
# CHANGE: Remove estimand == "ate" part
    if (estimand == "att") {
        if (!is.null(start)) {
            mean(Y[A == 1]) - rev(extrema.md(A, X, Y, gamma, "missing", reg.adjust, - start))
        } else {
            mean(Y[A == 1]) - rev(extrema.md(A, X, Y, gamma, "missing", reg.adjust))
        }
    }

}
