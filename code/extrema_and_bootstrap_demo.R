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

# Dan's code for plotting amplification
plotAmplificationBeta <- function(Z, X, Y, bound, Lambda, w, num_cov = 8, num_label = 2, in_data){
  ####################################################################
  # Plot amplification of bias = imbalance in U * \beta_u
  #
  # Standardize observed covariates for control units
  # Compute max regression coefficient among standardized observed covariates for control units
  #
  # Standardize observed covariates for all units
  # Compute maximum imbalance for standardized observed covariates before and after weighting
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   bound: if att >= 0, bound = upper bound on bias = sup h \hat{mu}_0^{(h)} - \hat{mu}_0
  #          if att < 0, bound = abs value of lower bound on bias = inf h \hat{mu}_0^{(h)} - \hat{mu}_0
  #   Lambda: sensitivity parameter
  #   w: estimated weights
  #   num_cov: number of observed covariates with highest values of \beta_u to plot. Default = 8
  #   num_label: number of observed covariates to include on plot. covariates with top values of 
  #              pre-weighting imbal * beta_u
  #
  # output:
  #   [[1]]: plot: contour plot of amplification of bias = imbalance in U * \beta_u
  #   [[2]]: df with regression coefficients and imbalance for each transformed
  #          observed covariate
  #
  ####################################################################
  
  ####################################################
  # Standardize observed covariates for control units
  ####################################################
  
  # get X for control units
  X_control <- X[Z == 0,]
  
  # standardize X for control units
  # center and make var = sd = 1
  X_control_stnd <- apply(X_control, MARGIN = 2, FUN = function(x) {(x- mean(x))/sd(x)})
  
  #################################################################
  # Compute max coefficient among standardized observed covariates
  #################################################################
  
  # max of coefficients for observed except intercept
  # each coefficient is with one variable standardized (U) and rest of variables not standardized
  #coeffs <- numeric()
  #for (var in 1:ncol(X_control)) {
  #  coeffs[var] <- lm(Y[Z==0]~X_control[,-var] + X_control_stnd[, var])$coef["X_control_stnd[, var]"]
  #}
  
  # equivalent to above is running on all standardized vars at same time
  coeffs <- lm(Y[Z==0]~X_control_stnd)$coef[-1]
  max_betau_01 <- max(abs(coeffs))
  
  ########################################################
  # Compute maximum imbalance for standardized covariates
  ########################################################
  
  X_stnd <- apply(X, MARGIN = 2, FUN = function(x) {(x - mean(x))/sd(x)})
  imbal_stnd <- colMeans(X_stnd[Z == 1,]) - colMeans(X_stnd[Z == 0,])
  max_imbal_stnd <- max(abs(imbal_stnd))
  
  ####################################################
  # Compute imbalance in covariates after weighting
  ####################################################
  
  #imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(Z))
  #imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(1-Z))
  imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(w))
  
  ####################################################
  # Get coordinates for strongest observed covariates to plot
  ####################################################
  # get coefficients
  coeff_df <- data.frame(
    covar = sub("X_control_stnd", "", names(coeffs)),
    coeff = abs(as.numeric(coeffs)))
  # get imbal
  imbal_df <- data.frame(
    covar = names(imbal_stnd),
    imbal = abs(as.numeric(imbal_stnd)),
    imbal_wt = abs(as.numeric(imbal_stnd_weight)))
  # merge coefficients and imbalance, arrange from largest to smallest by imbal * beta_u
  strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>%
    dplyr::arrange(desc(coeff*imbal))
  
  #######
  # Plot
  #######
  
  # function to compute beta_u from imbalance for plot
  betauFun <- function(x) {
    bound/x
  }
  
  # place holder until change to approx bal
  #imbal_df$imbal_wt <- imbal_df$imbal * 0.4
  # strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>% dplyr::arrange(desc(coeff*imbal))
  #
  
  # Create region for observed covariates post-weighting 
  
  x_orig = strongest_cov_df[1:num_cov,]$imbal_wt
  y_orig = strongest_cov_df[1:num_cov,]$coeff
  
  x <- c(x_orig,0, 0, max(x_orig))
  y <- c(y_orig,0, max(y_orig),0)
  df_plot <- matrix(c(x,y), ncol = 2)
  
  hpts <- chull(df_plot)
  hpts <- c(hpts, hpts[1])
  
  X_cvx <- df_plot[hpts,]
  
  # Create region for observed covariates post-weighting
  
  # x_orig = strongest_cov_df[1:num_cov,]$imbal_wt
  # y_orig = strongest_cov_df[1:num_cov,]$coeff
  # 
  # # compute convex hull for points (imbalance post-weighting and beta_u)
  # X <- matrix(c(x_orig,y_orig), ncol = 2)
  # hpts <- chull(X)
  # hpts <- c(hpts, hpts[1])
  # 
  # # add points for x-axis (x val for min(y), 0), y-axis (0, y val for min(x)), and origin
  # X_cvx <-
  #   rbind(
  #     X[hpts, ][1:which.min(X[hpts,2 ]),],
  #     matrix(c(x_orig[which.min(y_orig)], 0, 0, 0, 0, y_orig[which.min(x_orig)]), ncol = 2),
  #     X[hpts, ][(which.min(X[hpts,2 ])+1):length(hpts),]
  #     )
  
  if (in_data == 'lalonde') {
    mycurve1 <- as.data.frame(curve(from=0.01207, to=3, betauFun)) 
  } else if (in_data == 'lalonde interactions') {
    mycurve1 <- as.data.frame(curve(from=0.051, to=4.55, betauFun)) 
  } else if (in_data == 'fish') {
    mycurve1 <- as.data.frame(curve(from=0.451, to=2, betauFun))
  }
  if (in_data %in% c('lalonde', 'lalonde interactions')) {
    #change lalonde labels
    strongest_cov_df$covar <- 
      with(strongest_cov_df,
           ifelse(covar == 're75', "1975 real earnings",
                  ifelse(covar == 're74', "1974 real earnings", covar)))
    
    plot <-
      ggplot() +
      geom_line(data=mycurve1,aes(x=x,y=y, colour = "error")) +
      #stat_function(fun = betauFun, aes(colour = "error")) +
      geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal, y = coeff, colour = "unweighted observed covs")) +
      geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal_wt, y = coeff, colour = "weighted observed covs")) +
      geom_hline(aes(yintercept=max(strongest_cov_df[1:num_cov,"coeff"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
      geom_vline(aes(xintercept=max(strongest_cov_df[1:num_cov,"imbal"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
      geom_path(aes(x = X_cvx[,1], y = X_cvx[,2], colour ="weighted observed covs"), alpha = 0.5) + 
      geom_polygon(aes(x = X_cvx[,1], y = X_cvx[,2]), fill = "red", alpha = 0.2) +
      geom_text(data=strongest_cov_df[1:num_label,], aes(x = imbal, y = coeff, colour = "unweighted observed covs",label=covar),
                size = 3, hjust = 0.1, vjust = 2.2) +
      #annotate("rect", xmin = 0, xmax = max(strongest_cov_df[1:num_cov,]$imbal_wt), 
      #          ymin = 0, ymax = max(strongest_cov_df[1:num_cov,]$coeff), alpha = .3) +
      scale_x_continuous(name = "absolute standardized imbalance in U", limits = c(0, max(max_imbal_stnd,3))) +
      scale_y_continuous(name = TeX("absolute $\\beta_u$"), limits = c(0, max(max_betau_01, betauFun(0.21)))) +
      ggtitle("B") +
      #ggtitle(TeX(paste0("$\\beta_u$ vs. imbalance for  $\\Lambda$ = ", Lambda,", ", in_data, " data"))) +
      scale_colour_manual(values = c("black","#00BFC4" , "#F8766D"), 
                          breaks = c("error", "unweighted observed covs", "weighted observed covs")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = "none",
            axis.title.x = element_text(size = 10))
    #theme(legend.title = element_blank(), legend.position = "bottom")
  } else if (in_data == 'fish') {
    
    #change fish labels
    strongest_cov_df$covar <- 
      with(strongest_cov_df,
           ifelse(covar == 'race6', "non-Hispanic Asian",covar))
    
    plot <-
      ggplot() +
      geom_line(data=mycurve1,aes(x=x,y=y, colour = "error")) +
      #stat_function(fun = betauFun, aes(colour = "error")) +
      geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal, y = coeff, colour = "unweighted observed covs")) +
      geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal_wt, y = coeff, colour = "weighted observed covs")) +
      geom_hline(aes(yintercept=max(strongest_cov_df[1:num_cov,"coeff"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
      geom_vline(aes(xintercept=max(strongest_cov_df[1:num_cov,"imbal"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
      geom_path(aes(x = X_cvx[,1], y = X_cvx[,2], colour ="weighted observed covs"), alpha = 0.5) + 
      geom_polygon(aes(x = X_cvx[,1], y = X_cvx[,2]), fill = "red", alpha = 0.2) +
      geom_text(data=strongest_cov_df[1,], aes(x = imbal, y = coeff, colour = "unweighted observed covs",label=covar),
                size = 3, position = position_nudge(x = 0.48, y = 0.12)) +
      geom_text(data=strongest_cov_df[num_label,], aes(x = imbal, y = coeff, colour = "unweighted observed covs",label=covar),
                size = 3, position = position_nudge(x = 0.235, y = 0)) +
      #annotate("rect", xmin = 0, xmax = max(strongest_cov_df[1:num_cov,]$imbal_wt), 
      #          ymin = 0, ymax = max(strongest_cov_df[1:num_cov,]$coeff), alpha = .3) +
      scale_x_continuous(name = "absolute standardized imbalance in U", limits = c(0, max(max_imbal_stnd,2))) +
      scale_y_continuous(name = TeX("absolute $\\beta_u$"), limits = c(0, max(max_betau_01, 4))) +
      ggtitle("B") +
      #ggtitle(TeX(paste0("$\\beta_u$ vs. imbalance for  $\\Lambda$ = ", Lambda,", ", in_data, " data"))) +
      scale_colour_manual(values = c("black","#00BFC4" , "#F8766D"), 
                          breaks = c("error", "unweighted observed covs", "weighted observed covs")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = "none")
    #theme(legend.title = element_blank(), legend.position = "bottom")
  }
  
  return(list(plot, strongest_cov_df %>% dplyr::mutate_if(is.numeric, round, digits = 3)))
}



