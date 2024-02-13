# Zhao's code adapted for disparity decomp


getExtrema <- function(G, Y, gamma = 0, w, estimand = "point") {
  estimand <- match.arg(estimand, c("point", "reduction", "residual"))
  # assume observed rmpw weights already computed
  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])
  mu_10 <- sum(w * Y * G) / sum(w * G)
  # message("mu1: ", round(mu1, 2))
  # message("mu0: ", round(mu0, 2))
  
  state <- switch(estimand, point = round(mu_10, 3), 
                  reduction = round(mu1 - mu_10, 3), 
                  residual = round(mu_10 - mu0, 3))
  
  message("Parameter of interest: ", estimand)
  message("Estimate: ", state)
  
  w <- w[G != 0]
  Y <- Y[G != 0]
  G <- G[G != 0]
  
  w <- w[order(-Y)]
  Y <- Y[order(-Y)]
  G <- G[order(-Y)]
  
  ## maximization
  num.each.low <- G * Y * (exp(-gamma) * w)
  num.each.up <- G * Y * (exp(gamma) * w)
  num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)
  
  den.each.low <- G * (exp(-gamma) * w)
  den.each.up <- G * (exp(gamma) * w)
  den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)
  
  maximum <- max(num / den)
  ## print(den[which.max(num/den)] / n)
  
  ## minimization
  num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
  minimum <- min(num / den)
  ## print(den[which.min(num/den)] / n)
  
  out <- c(minimum, maximum)
  switch(estimand, 
         point = out,
         reduction = mu1 - out,
         residual = out - mu0)
}


boostrapCI <- function(G, Y, gamma = 0, w, alpha = 0.05, estimand = "point", 
                        parallel = TRUE, B = 1000, warm.start = FALSE) {
  
  estimand <- match.arg(estimand, c("point", "reduction", "residual"))
  # assume observed rmpw weights already computed
  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])
  mu_10 <- sum(w * Y * G) / sum(w * G)
  # message("mu1: ", round(mu1, 2))
  # message("mu0: ", round(mu0, 2))
  
  state <- switch(estimand, point = round(mu_10, 3), 
                  reduction = round(mu1 - mu_10, 3), 
                  residual = round(mu_10 - mu0, 3))
  
  message("Parameter of interest: ", estimand)
  message("Estimate: ", state)
  
  no.cores <- max(1, ifelse(parallel, detectCores(), 1))
  n <- length(G)
  
  # if (warm.start) {
  #   start <- glm(A ~ X, family = "binomial")$coefs
  # } else {
  #   start <- NULL
  # }
  
  out <- mclapply(1:B, function(iter) {
    s <- sample(1:n, n, TRUE);
    res <- tryCatch(getExtrema(G[s], Y[s], gamma, w[s], estimand),
                    error = function(e) {print(e)});
    res
  }, mc.cores = no.cores)
  
  out <- do.call(rbind, out)
  
  c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))
}

