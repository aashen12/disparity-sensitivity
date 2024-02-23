
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
  
  #subset_lalonde1 <- sample(1:185, size = 185, replace = F)
  #subset_lalonde2 <- sample(186:length(Z), size = 500, replace = F)
  #subset_lalonde <- c(subset_lalonde1, subset_lalonde2)
  
  #Z <- Z[subset_lalonde]
  #Y <- Y[subset_lalonde]
  #X <- X[subset_lalonde,]
  
  # Estimate weights
  #weights_x <- ebalance(Treatment = Z, X = X)$w
  #weights_x <- balancer::multilevel_qp(X = X, trt = Z, Z = rep(1, nrow(X)), exact_global = F)$weights[Z==0]
  data_frame = as.data.frame(cbind(Z, X, Y))
  # Define treatment indicator and
  z_ind = "Z"
  # moment covariates
  bal = list()
  bal$bal_cov = colnames(X)
  # Set tolerances
  bal$bal_tol = 0.02
  bal$bal_std = "group"
  # Solve for the Average Treatment Effect on the Treated, ATT (default)
  bal$bal_alg = FALSE
  system.time(
  sbwatt_object <- sbw::sbw(dat = data_frame, ind = z_ind, out = "Y", bal = bal)
  )
  weights_x <- sbwatt_object$dat_weights[sbwatt_object$dat_weights$Z == 0,]$sbw_weights
  
  # Compute \mu_0 hat
  Y_0_for_control <- Y[Z == 0]
  mu_0_hat <- sum(Y_0_for_control * weights_x)
  
  # compute bounds
  bounds <- rev(extrema.md(Z, X, Y, gamma = log(Lambda), "missing", reg.adjust = FALSE)) - mu_0_hat
  
  # return bounds and mu_0_hat
  out <- list(bounds, mu_0_hat, weights_x)
  
  return(out)
}
