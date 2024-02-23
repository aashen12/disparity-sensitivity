
getBiasBounds <- function(Z, X, Y, Lambda) {
  ####################################################################
  # Get bounds on bias = \mu^*_0 - \hat{mu}_0
  # lower bound = inf h \hat{mu}_0^{(h)} - \hat{mu}_0
  # upper bound = inf h \hat{mu}_0^{(h)} - \hat{mu}_0
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   Lambda: sensitivity parameter
  #
  # output:
  #   bounds: out[[1]]
  #     bounds[1]: upper bound = sup h \hat{mu}_0^{(h)} - \hat{mu}_0
  #     bounds[2]: lower bound = inf h \hat{mu}_0^{(h)} - \hat{mu}_0
  #   mu_0_hat: out[[2]]. estimate of E[Y(0)|Z==1]
  #   eb.out$w: out[[3]]. Estimated balancing weights
  #
  ####################################################################
  
  # Estimate weights
  #weights_x <- ebalance(Treatment = Z, X = X)$w
  #weights_x <- balancer::multilevel_qp(X = X, trt = Z, Z = rep(1, nrow(X)), exact_global = F)$weights[Z==0]
  tol <- 0.15
  sbw_out <- sbw_osqp(scale(X), Z, tol, verbose = TRUE)
  weights_x <- sbw_out$wts[Z==0]
  # Compute \mu_0 hat
  Y_0_for_control <- Y[Z == 0]
  #n_treat <- sum(Z)
  Y_0_for_control_times_w <- Y_0_for_control * weights_x
  #mu_0_hat <- (1/n_treat)*sum(Y_0_for_control_times_w)
  mu_0_hat <- (1/sum(weights_x))*sum(Y_0_for_control_times_w)
  
  # compute bounds
  bounds <- rev(extrema.md(Z, X, Y, gamma = log(Lambda), "missing", reg.adjust = FALSE)) - mu_0_hat
  
  # return bounds and mu_0_hat
  out <- list(bounds, mu_0_hat, weights_x)
  
  return(out)
}
