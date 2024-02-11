
#################################################
# Amplification: bias = imbalance in U * \beta_u
#################################################
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

plotCIs <- function(Z, X, Y, Lambda_vec, in_data){
  ####################################################################
  # Plot point estimate and confidence intervals for varying values of Lambda
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   Lambda_vec: values of sensitivity analysis parameter for which to plot results  
  #   input_data: name of data to plot
  #
  # output:
  #   [[1]]: plot: plot of sensitivity analysis results for varying values of Lambda
  #   [[2]]: table: table of sensitivity analysis results for varying values of Lambda
  #
  ####################################################################


#####
# balancing weights estimate:
#####

# point estimate
point1 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[1]))
point2 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[2]))
point3 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[3]))
point4 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[4]))
# conf int
CI1 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                   gamma = log(Lambda_vec[1]))
CI2 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                   gamma = log(Lambda_vec[2]))
CI3 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                   gamma = log(Lambda_vec[3]))
CI4 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                   gamma = log(Lambda_vec[4]))

boot_fig_df <- data.frame(Lambda = Lambda_vec,
                          point_min = c(point1[1],
                                        point2[1],
                                        point3[1],
                                        point4[1]),
                          point_max = c(point1[2],
                                        point2[2],
                                        point3[2],
                                        point4[2]),
                          conf_min = c(CI1[1],
                                       CI2[1],
                                       CI3[1],
                                       CI4[1]),
                          conf_max = c(CI1[2],
                                       CI2[2],
                                       CI3[2],
                                       CI4[2]),
                          grp = c("bal", 
                                  "bal",
                                  "bal", 
                                  "bal")
)

boot_fig_df$point_mean <- (boot_fig_df$point_min + boot_fig_df$point_max)/2

boot_fig_df$Lambda <- as.factor(boot_fig_df$Lambda)

plot <- 
  ggplot(boot_fig_df, aes(x = Lambda, y = point_mean)) +
  geom_point(size = 0) +
  geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.5) +
  geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.45) +
  #scale_x_continuous(breaks=c(1, 2, 5.33, 6)) +
  #ggtitle(paste0("sensitivity analysis results - ", in_data, " data")) +
  ylab("ATT") +
  xlab(TeX("$\\Lambda$")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# table

table_out <- data.frame(
  Lambda = boot_fig_df$Lambda,
  point_est = paste0("[",round(boot_fig_df$point_min,2), ", ", round(boot_fig_df$point_max,2),"]"),
  conf_int = paste0("[",round(boot_fig_df$conf_min,2), ", ", round(boot_fig_df$conf_max,2),"]")
)

colnames(table_out) <- c("Lambda", "point est", "95% conf int")

return(list(plot, table_out))
}

#################################################
# Amplification: bias = imbalance in U * \beta_u
#################################################
plotConceptualAmplification <- function(Z, X, Y, bound_low, bound_mid, bound_high, w, num_cov = 8, num_label = 0){
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
  #   scenario: "sensitive", "ambiguous", "robust"
  #   w: estimated weights
  #   num_cov: number of observed covariates with highest values of \beta_u to plot. Default = 8
  #   num_label: number of observed covariates to include on plot. covariates with top values of 
  #              pre-weighting imbal * beta_u
  #
  # output:
  #   plot: contour plot of amplification of bias = imbalance in U * \beta_u
  #
  ####################################################################
  
  ####################################################
  # Standardize observed covariates for control units
  ####################################################
  
  # get X for control units
  X_control <- X[Z == 0,]
  
  # standardize X for control units
  X_control_stnd <- apply(X_control, MARGIN = 2, FUN = function(x) {x/sd(x)})
  
  #################################################################
  # Compute max coefficient among standardized observed covariates
  #################################################################
  
  # max of coefficients for observed except intercept
  coeffs <- lm(Y[Z==0]~X_control_stnd)$coef[-1]
  max_betau_01 <- max(abs(coeffs))
  
  ########################################################
  # Compute maximum imbalance for standardized covariates
  ########################################################
  
  X_stnd <- apply(X, MARGIN = 2, FUN = function(x) {x/sd(x)})
  imbal_stnd <- colMeans(X_stnd[Z == 1,]) - colMeans(X_stnd[Z == 0,])
  max_imbal_stnd <- max(abs(imbal_stnd))
  
  ####################################################
  # Compute imbalance in covariates after weighting
  ####################################################
  
  #imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(Z))
  imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(1-Z))
  
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
  betauFun_low <- function(x) {
    bound_low/x
  }
  
  # function to compute beta_u from imbalance for plot
  betauFun_mid <- function(x) {
    bound_mid/x
  }
  
  # function to compute beta_u from imbalance for plot
  betauFun_high <- function(x) {
    bound_high/x
  }
  # place holder until change to approx bal

  imbal_df[imbal_df$covar == "income",]$imbal <- 1
  imbal_df[imbal_df$covar == "race1",]$imbal <- 0.5
  imbal_df[imbal_df$covar == "race6",]$imbal <- 0.6
  imbal_df[imbal_df$covar == "race3",]$imbal <- 0.24
  
  imbal_df$imbal_wt <- imbal_df$imbal * 0.4
  
  imbal_df[imbal_df$covar == "race6",]$imbal_wt <- 0.3
  imbal_df[imbal_df$covar == "income",]$imbal_wt <- 0.5
  
  strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>% dplyr::arrange(desc(coeff*imbal))
  
  
  strongest_cov_df[strongest_cov_df$covar == "race6",]$coeff <- 0.25
  strongest_cov_df[strongest_cov_df$covar == "race3",]$coeff <- 0.09
  
  # Create region for observed covariates post-weighting 
  
  x_orig = strongest_cov_df[1:num_cov,]$imbal_wt
  y_orig = strongest_cov_df[1:num_cov,]$coeff
  
  x <- c(x_orig,0, 0, max(x_orig))
  y <- c(y_orig,0, max(y_orig),0)
  #x <- c(x_orig,0, 0+0.0001, x_orig[which.min(y_orig)]-0.0001)
  #y <- c(y_orig,0,y_orig[which.min(x_orig)]-0.0001,0 +0.0001)
  #y <- y[order(x)]
  #x <- x[order(x)]
  df_plot <- matrix(c(x,y), ncol = 2)
  
  # compute convex hull for points (imbalance post-weighting and beta_u)
  #X <- matrix(c(x_orig,y_orig), ncol = 2)
  hpts <- chull(df_plot)
  hpts <- c(hpts, hpts[1])
  
  # add points for x-axis (x val for min(y), 0), y-axis (0, y val for min(x)), and origin
  #X_cvx <-
  #  rbind(
  #    X[hpts, ][1:which.min(X[hpts,2 ]),],
  #    matrix(c(x_orig[which.min(y_orig)], 0, 0, 0, 0, y_orig[which.min(x_orig)]), ncol = 2),
  #    X[hpts, ][(which.min(X[hpts,2 ])+1):length(hpts),]
  #  )
  X_cvx <- df_plot[hpts,]
  # greatest convex minorant
  #gg = gcmlcm(x,y)
  # least concave majorant
  #ll = gcmlcm(x,y, type="lcm")
  
  #x.knots = c(gg$x.knots,rev(ll$x.knots)) 
  #y.knots = c(gg$y.knots,rev(ll$y.knots))
  #lines(x.knots, y.knots, col=5, lwd=2)
  
  plot <-
    ggplot() +
    stat_function(fun = betauFun_low, aes(colour = "bound a")) +
    stat_function(fun = betauFun_mid, aes(colour = "bound b")) +
    stat_function(fun = betauFun_high, aes(colour = "bound c")) +
    geom_text(aes(label = "sensitive to observed", colour = "bound a", x = 1.1, y = 0), hjust = -.1, size = 3) +
    geom_text(aes(label = "sensitive to unobserved", colour = "bound b", x = 1.22, y = 0.12), hjust = -.1, size = 3) +
    geom_text(aes(label = "robust", colour = "bound c", x = 1.5, y = 0.28), hjust = -.1, size = 3) +
    geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal, y = coeff, colour = "unweighted observed covs")) +
    geom_point(data = strongest_cov_df[1:num_cov,], aes(x = imbal_wt, y = coeff, colour = "weighted observed covs")) +
    geom_hline(aes(yintercept=max(strongest_cov_df[1:num_cov,"coeff"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
    geom_vline(aes(xintercept=max(strongest_cov_df[1:num_cov,"imbal"]), colour = "unweighted observed covs"),linetype="dashed", alpha = 0.4) +
    geom_path(aes(x = X_cvx[,1], y = X_cvx[,2], colour ="weighted observed covs"), alpha = 0.5) + 
    geom_polygon(aes(x = X_cvx[,1], y = X_cvx[,2]), fill = "red", alpha = 0.2) +
    #geom_text(data=strongest_cov_df[1:num_label,], aes(x = imbal, y = coeff, colour = "observed covs",label=covar),
              #size = 3, position = position_nudge(x = 0.3)) +
    #annotate("rect", xmin = 0, xmax = max(strongest_cov_df[1:num_cov,]$imbal_wt), 
    #          ymin = 0, ymax = max(strongest_cov_df[1:num_cov,]$coeff), alpha = .3) +
    scale_x_continuous(name = "absolute standardized imbalance in U", limits = c(0, max(max_imbal_stnd,2))) +
    scale_y_continuous(name = TeX("absolute $\\beta_u$"), limits = c(0, 0.75)) +
    #ggtitle(TeX(paste0("$\\beta_u$ vs. imbalance, scenario: ", scenario))) +
    scale_colour_manual(values = c("black","#00BA38","#C77CFF","#00BFC4" , "#F8766D"), 
                        breaks = c("", "", "", "unweighted observed covs", "weighted observed covs"),
                        labels = c("", "", "", "unweighted observed covs", "weighted observed covs")) +
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "bottom")
  
  return(plot)
}

plotCIsAnim <- function(Z, X, Y, Lambda_vec, in_data){
  ####################################################################
  # Plot point estimate and confidence intervals for varying values of Lambda
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   Lambda_vec: values of sensitivity analysis parameter for which to plot results  
  #   input_data: name of data to plot
  #
  # output:
  #   [[1]]: plot: plot of sensitivity analysis results for varying values of Lambda
  #   [[2]]: table: table of sensitivity analysis results for varying values of Lambda
  #
  ####################################################################
  
  
  #####
  # balancing weights estimate:
  #####
  
  # point estimate
  point1 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[1]))
  point2 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[2]))
  point3 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[3]))
  point4 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[4]))
  # conf int
  CI1 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[1]))
  CI2 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[2]))
  CI3 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[3]))
  CI4 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[4]))
  
  boot_fig_df <- data.frame(Lambda = Lambda_vec,
                            point_min = c(point1[1],
                                          point2[1],
                                          point3[1],
                                          point4[1]),
                            point_max = c(point1[2],
                                          point2[2],
                                          point3[2],
                                          point4[2]),
                            conf_min = c(CI1[1],
                                         CI2[1],
                                         CI3[1],
                                         CI4[1]),
                            conf_max = c(CI1[2],
                                         CI2[2],
                                         CI3[2],
                                         CI4[2]),
                            grp = c("bal", 
                                    "bal",
                                    "bal", 
                                    "bal")
  )
  
  boot_fig_df$point_mean <- (boot_fig_df$point_min + boot_fig_df$point_max)/2
  
  boot_fig_df$Lambda <- as.factor(boot_fig_df$Lambda)
  
  alpha_list <- c(1, 0, 0, 0)
  color_rule <- ifelse(alpha_list == 0, NA, "black")

  plot1 <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean)) +
    geom_point(size = 0, alpha = alpha_list) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.5, color = color_rule) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.5, color = color_rule) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = 'none')
  
  alpha_list <- c(1, 1, 0, 0)
  color_rule <- ifelse(alpha_list == 0, NA, "black")
  
  plot2 <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean)) +
    geom_point(size = 0, alpha = alpha_list) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.5, color = color_rule) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.5, color = color_rule) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = 'none')
  
  alpha_list <- c(1, 1, 1, 0)
  color_rule <- ifelse(alpha_list == 0, NA, "black")
  
  plot3 <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean)) +
    geom_point(size = 0, alpha = alpha_list) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.5, color = color_rule) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.5, color = color_rule) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = 'none')
  
  
  alpha_list <- c(1, 1, 1, 1)
  color_rule <- ifelse(alpha_list == 0, NA, "black")
  
  plot4 <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean)) +
    geom_point(size = 0, alpha = 0) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.5) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.5) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = 'none')
  
  # table
  
  table_out <- data.frame(
    Lambda = boot_fig_df$Lambda,
    point_est = paste0("[",round(boot_fig_df$point_min,2), ", ", round(boot_fig_df$point_max,2),"]"),
    conf_int = paste0("[",round(boot_fig_df$conf_min,2), ", ", round(boot_fig_df$conf_max,2),"]")
  )
  
  colnames(table_out) <- c("Lambda", "point est", "95% conf int")
  
  return(list(plot1,plot2, plot3, plot4))
}

plotCICompare <- function(Z, X, Y, Lambda_vec, in_data){
  ####################################################################
  # Plot point estimate and confidence intervals for varying values of Lambda
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   Lambda_vec: values of sensitivity analysis parameter for which to plot results  
  #   input_data: name of data to plot
  #
  # output:
  #   [[1]]: plot: plot of sensitivity analysis results for varying values of Lambda
  #   [[2]]: table: table of sensitivity analysis results for varying values of Lambda
  #
  ####################################################################
  
  
  #####
  # balancing weights estimate:
  #####
  Lambda_vec <- c(1, 1.65, 2.72, 4, 7.39)
  source("../r_functions/DG_qreg_ATT.R")
  library(quantreg)
  
  # point estimate
  point1 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[1]))
  point2 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[2]))
  point3 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[3]))
  point4 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[4]))
  point5 <- extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[5]))
  # conf int
  CI1 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[1]))
  CI2 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[2]))
  CI3 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[3]))
  CI4 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[4]))
  CI5 <- bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[5]))
  
  # Zhao
  # point estimate
  point1_zsb <- bootsens::extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[1]))
  point2_zsb <- bootsens::extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[2]))
  point3_zsb <- bootsens::extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[3]))
  point4_zsb <- bootsens::extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                       gamma = log(Lambda_vec[4]))
  point5_zsb <- bootsens::extrema.os(Z, X, Y, estimand = "att", reg.adjust = FALSE, 
                                     gamma = log(Lambda_vec[5]))
  # conf int
  CI1_zsb <- bootsens::bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[1]))
  CI2_zsb <- bootsens::bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[2]))
  CI3_zsb <- bootsens::bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[3]))
  CI4_zsb <- bootsens::bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                     gamma = log(Lambda_vec[4]))
  CI5_zsb <- bootsens::bootsens.os(Z, X, Y, alpha = 0.05, estimand = "att", reg.adjust = FALSE, 
                                   gamma = log(Lambda_vec[5]))
  
  # D&G
  # point estimate
  point1_dg <- run_bootstrap_iter(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[1], bootstrap = FALSE)
  point2_dg <- run_bootstrap_iter(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[2], bootstrap = FALSE)
  point3_dg <- run_bootstrap_iter(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[3], bootstrap = FALSE)
  point4_dg <- run_bootstrap_iter(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[4], bootstrap = FALSE)
  
  # conf int
  CI1_dg <- obtainCIsDG(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[1], nboot = 1000)
  CI2_dg <- obtainCIsDG(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[2], nboot = 1000)
  CI3_dg <- obtainCIsDG(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[3], nboot = 1070)
  CI4_dg <- obtainCIsDG(df_sample = data.frame(Y,Z,X), lambda = Lambda_vec[4], nboot = 1350)
  
  boot_fig_df <- data.frame(Lambda = c(Lambda_vec,Lambda_vec,Lambda_vec[1:4]),
                            point_min = c(point1[1],
                                          point2[1],
                                          point3[1],
                                          point4[1],
                                          point5[1],
                                          point1_zsb[1],
                                          point2_zsb[1],
                                          point3_zsb[1],
                                          point4_zsb[1],
                                          point5_zsb[1],
                                          point1_dg[1],
                                          point2_dg[1],
                                          point3_dg[1],
                                          point4_dg[1]),
                            point_max = c(point1[2],
                                          point2[2],
                                          point3[2],
                                          point4[2],
                                          point5[2],
                                          point1_zsb[2],
                                          point2_zsb[2],
                                          point3_zsb[2],
                                          point4_zsb[2],
                                          point5_zsb[2],
                                          point1_dg[2],
                                          point2_dg[2],
                                          point3_dg[2],
                                          point4_dg[2]),
                            conf_min = c(CI1[1],
                                         CI2[1],
                                         CI3[1],
                                         CI4[1],
                                         CI5[1],
                                         CI1_zsb[1],
                                         CI2_zsb[1],
                                         CI3_zsb[1],
                                         CI4_zsb[1],
                                         CI5_zsb[1],
                                         CI1_dg[1],
                                         CI2_dg[1],
                                         CI3_dg[1],
                                         CI4_dg[1]),
                            conf_max = c(CI1[2],
                                         CI2[2],
                                         CI3[2],
                                         CI4[2],
                                         CI5[2],
                                         CI1_zsb[2],
                                         CI2_zsb[2],
                                         CI3_zsb[2],
                                         CI4_zsb[2],
                                         CI5_zsb[2],
                                         CI1_dg[2],
                                         CI2_dg[2],
                                         CI3_dg[2],
                                         CI4_dg[2]),
                            grp = rep(c("bal","zsb","dg"), times = c(rep(5,2),4))
  )
  
  boot_fig_df$Lambda <- 
    with(boot_fig_df,
         ifelse(grp == 'zsb', Lambda - 0.22,
                ifelse(grp == 'dg', Lambda + 0.22, Lambda
                )))
  
  #boot_fig_df[boot_fig_df$Lambda == '7.19',]$Lambda <- 7.29
  #boot_fig_df[boot_fig_df$Lambda == '7.39',]$Lambda <- 7.49
    
  boot_fig_df$point_mean <- (boot_fig_df$point_min + boot_fig_df$point_max)/2
  
  #boot_fig_df$Lambda <- as.factor(boot_fig_df$Lambda)
  
  boot_fig_df$grp <- factor(boot_fig_df$grp, 
                            levels = c("zsb", "bal", "dg"))
  
  plot <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean, color = grp)) +
    geom_point(size = 0) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.2) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.2) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    scale_x_continuous(breaks=Lambda_vec) +
    #ggtitle(paste0("sensitivity analysis comparison - ", in_data, " data")) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.position = "none")
    #theme(legend.title = element_blank())
  
  
  
  boot_fig_df <- data.frame(Lambda = Lambda_vec,
                            point_min = c(point1[1],
                                          point2[1],
                                          point3[1],
                                          point4[1],
                                          point1_zsb[1],
                                          point2_zsb[1],
                                          point3_zsb[1],
                                          point4_zsb[1]),
                            point_max = c(point1[2],
                                          point2[2],
                                          point3[2],
                                          point4[2],
                                          point1_zsb[2],
                                          point2_zsb[2],
                                          point3_zsb[2],
                                          point4_zsb[2]),
                            conf_min = c(CI1[1],
                                         CI2[1],
                                         CI3[1],
                                         CI4[1],
                                         CI1_zsb[1],
                                         CI2_zsb[1],
                                         CI3_zsb[1],
                                         CI4_zsb[1]),
                            conf_max = c(CI1[2],
                                         CI2[2],
                                         CI3[2],
                                         CI4[2],
                                         CI1_zsb[2],
                                         CI2_zsb[2],
                                         CI3_zsb[2],
                                         CI4_zsb[2]),
                            grp = rep(c("bal","zsb"), times = rep(4,2))
  )
  
  boot_fig_df$Lambda <- 
    with(boot_fig_df,
         ifelse(grp == 'zsb', Lambda - 0.15,
                ifelse(grp == 'bal', Lambda + 0.15, Lambda
                )))
  
  boot_fig_df$point_mean <- (boot_fig_df$point_min + boot_fig_df$point_max)/2
  
  #boot_fig_df$Lambda <- as.factor(boot_fig_df$Lambda)
  
  boot_fig_df$grp <- factor(boot_fig_df$grp, 
                            levels = c("zsb", "bal"))
  
  plot <- 
    ggplot(boot_fig_df, aes(x = Lambda, y = point_mean, color = grp)) +
    geom_point(size = 0) +
    geom_errorbar(aes(ymin = point_min, ymax = point_max), width = 0.18) +
    geom_errorbar(aes(ymin = conf_min, ymax = conf_max), linetype = "dashed", width = 0.18) +
    geom_hline(yintercept = 0, alpha = 0.45) +
    scale_x_continuous(breaks=Lambda_vec) +
    #ggtitle(paste0("sensitivity analysis comparison - ", in_data, " data")) +
    ylab("ATT") +
    xlab(TeX("$\\Lambda$")) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  
  
  
  
  # table
  
  table_out <- data.frame(
    Lambda = boot_fig_df$Lambda,
    point_est = paste0("[",round(boot_fig_df$point_min,2), ", ", round(boot_fig_df$point_max,2),"]"),
    conf_int = paste0("[",round(boot_fig_df$conf_min,2), ", ", round(boot_fig_df$conf_max,2),"]")
  )
  
  colnames(table_out) <- c("Lambda", "point est", "95% conf int")
  
  return(list(plot, table_out))
}
