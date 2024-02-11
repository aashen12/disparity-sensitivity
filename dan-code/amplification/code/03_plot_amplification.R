#' ---
#' title: "Sensitivity Analysis for Balancing Weights"
#' output: html_document
#' ---

#' # Amplification

#+ load_libs_data, include = FALSE
rm(list = ls())

library(foreign)
library(tidyverse)
#library(ebal)
library(parallel)
library(scales)
library(kableExtra)
library(latex2exp)
library(ppcor)
library(fdrtool)
library(balancer)
#library(sbw)

# load input df
input_df <- readRDS("../inputs/input_df.rds")

# load balancing weights sensitivity functions
if (input_df$data == 'lalonde interactions') {
  source("../r_functions/new_ipw_BW_interactions_tol.R")
  source("../r_functions/amplification_functions_interactions_tol.R")
} else {
  source("../r_functions/new_ipw_BW.R")
  source("../r_functions/amplification_functions.R")
}
source("../r_functions/tmp_sbw.R")

# load amplification functions
source("../r_functions/eli_sbw_plot_functions_new.R")

######################
# load data
######################

# treatment vector
Z <- readRDS(file = paste0("../", input_df$data, "/intermediate/Z_", input_df$data, ".rds"))
# design matrix
X <- readRDS(file = paste0("../", input_df$data, "/intermediate/X_", input_df$data, ".rds"))
# outcome vector
Y <- readRDS(file = paste0("../", input_df$data, "/intermediate/Y_", input_df$data, ".rds"))

#+ compute_bounds, include = FALSE
##################
# Compute lambda*
##################
#extrema.os(Z,X,Y, estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.08), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.18), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.03), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.03), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.02), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.02), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.01), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.01), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.01), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.05), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.04), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.03), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(1.02), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(5.9), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(5.4), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(5.45), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log(5.5), estimand = 'att')
#bootsens.os(Z,X,Y, gamma = log([INSERT]), estimand = 'att')
#############
# Set Lambda
#############
if (input_df$data == 'rhc') {
  Lambda <- 1.08
} else if (input_df$data == 'nhanes') {
  Lambda <- 1.18
} else if (input_df$data == 'lalonde') {
  Lambda <- 1.01
 # Lambda <- 1.02
} else if (input_df$data == 'lalonde interactions') {
  Lambda <- 1.05
} else if (input_df$data == 'school') {
  
} else if (input_df$data == 'fish') {
  #Lambda <- 5.9
  Lambda <- 5.5
}

####################################
# Set Lambda range for results plot
####################################
if (input_df$data == 'rhc') {
  Lambda <- 1.08
} else if (input_df$data == 'nhanes') {
  Lambda <- 1.18
} else if (input_df$data == 'lalonde') {
  Lambda_vec <- c(1, 1.01, 1.25, 1.5)
  # Lambda_vec <- c(1, 1.02, 1.25, 1.5)
} else if (input_df$data == 'lalonde interactions') {
  Lambda_vec <- c(1, 1.05, 1.25, 1.5)
} else if (input_df$data == 'school') {
  
} else if (input_df$data == 'fish') {
  Lambda_vec <- c(1, 3, 5.5, 7)
}

######################
# compute bias bounds
######################
out_bounds <- getBiasBounds(Z, X, Y, Lambda)
bounds <- out_bounds[[1]]
# upper and lower bounds
upper <- bounds[1]
lower <- bounds[2]

# mu_0_hat: estimate of E[Y(0)|Z==1]
mu_0_hat <- out_bounds[[2]]

# weights
w <- out_bounds[[3]]

# Estimate ATT:
att_est <- mean(Y[Z == 1]) - mu_0_hat

# Use lower or upper bound depending on if ATT estimate is positive or negative
#if (att_est >= 0) {
#  bound <- upper
#} else if (att_est < 0) {
#  bound <- abs(lower)
#}
# Use max abs value of bounds
bound <- max(abs(upper), abs(lower))

#' Amplification of bias = imbalance in $U \times \beta_u$

#+ ampl_beta, include = TRUE, echo = FALSE, warning = FALSE, cache = TRUE
#################################################
# Amplification: bias = imbalance in U * \beta_u
#################################################

amp_beta <- plotAmplificationBeta(Z, X, Y, bound, Lambda, w, num_cov = ncol(X), num_label = 2, in_data = input_df$data)

# plot
amp_beta[[1]]

#' \begin{itemize}
#' \item \textbf{bound:} 
#' \begin{itemize}
#' \item if estimated ATT is positive, bound = $\Big(\underset{h\in \mathcal{H}(\Lambda)}{\sup}\hat{\mu}^{(h)}_0 \Big) - \hat{\mu}_0$
#' \item if estimated ATT is negative, bound = $\Big(\underset{h\in \mathcal{H}(\Lambda)}{\inf}\hat{\mu}^{(h)}_0 \Big) - \hat{\mu}_0$
#' \end{itemize}
#' \item We consider $U \in [0,1]$, so we transform each observed covariate as follows:
#' \begin{itemize}
#' \item Make min = 0: subtract min value of covariate
#' \item Make max = 1: divide by max of shifted covariate
#' \end{itemize}
#' \item \textbf{max $\beta$ obs:} max absolute value of coefficients of transformed covariates from OLS of Y on transformed covariates for control units.
#' \item \textbf{max imbal obs:} max absolute value of difference in means of transformed covariates before weighting between treatment and control.
#' \item \textbf{top $\beta$ obs:} coefficient and imbalance for specified number of observed covariates sorted by descending coefficient value
#' \end{itemize}

#+ ampl_beta_table, include = TRUE, echo = FALSE, warning = FALSE
# table

if (input_df$data == 'fish') {
  kable(amp_beta[[2]] %>% rename(covariate = covar, coefficient = coeff, imbalance = imbal, `post-weighting imbalance` = imbal_wt) %>% 
          dplyr::mutate_if(is.numeric, round, digits = 3)) %>%
    kable_styling(position = "center")
} else if (input_df$data == 'lalonde'){
  kable(amp_beta[[2]] %>% rename(covariate = covar, coefficient = coeff, imbalance = imbal, `post-weighting imbalance` = imbal_wt) %>% 
          dplyr::mutate_if(is.numeric, round, digits = 2)) %>%
    kable_styling(position = "center")
}

#' Confidence and point estimate intervals

#+ conf_int, include = TRUE, echo = FALSE, warning = FALSE

conf_int <- plotCIs(Z,X,Y, Lambda_vec, in_data = input_df$data)

# plot
conf_int[[1]]

#conf_int_anim <- plotCIsAnim(Z,X,Y, Lambda_vec, in_data = input_df$data)
#conf_int_anim[[1]]
#conf_int_anim[[2]]
#conf_int_anim[[3]]
#conf_int_anim[[4]]

# table
kable(conf_int[[2]] %>% rename(`point estimate` = `point est`, `95% confidence interval` = `95% conf int`)) %>%
  kable_styling(position = "center")

#' Comparison with Zhao et al. and Dorn & Guo
source("../r_functions/DG_qreg_ATT.R")

#+ conf_int_zhao_dg, include = TRUE, echo = FALSE, warning = FALSE

conf_int <- plotCIs(Z,X,Y, Lambda_vec, in_data = input_df$data)
plotCICompare

# plot
conf_int[[1]]
