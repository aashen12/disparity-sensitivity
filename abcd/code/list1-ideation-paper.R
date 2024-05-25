# Performs amplification of the sensitivity analysis
# Analysis of list 1
rm(list = ls())
# devtools::install_github("https://github.com/aashen12/decompsens")
# devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))


library(tidyverse)
library(decompsens)
library(parallel)
library(doParallel)
library(boot)
library(jointVIP)
library(latex2exp)

options(na.action='na.pass')
source("functions.R")
set.seed(122357)

numCores <- parallel::detectCores()

doParallel::registerDoParallel(numCores)

Z_method <- "worry_upset"
# "better_worry",
# "smile",
# "better_upset",
# "love",
# "easy_talk
# aggregate
# worry_upset

outcome <- "ideation" # ideation or attempt

df_x <- read_csv(paste0("../data/list1_X_", Z_method, "_", outcome, ".csv"))
df_yz <- read_csv(paste0("../data/list1_YZGW_", Z_method, "_", outcome, ".csv"))

mediators <- c("src_subject_id")
mediators <- c("src_subject_id", "family_mental_health")
# allowable_covs <- c("age", "sex", "sib_order", "sib_num")
allowable_covs <- c("age", "sex")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df <- df_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

df <- df_non_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

interact <- FALSE




G <- df_yz$sex_min
Z <- df_yz$parent_accept
Y <- df_yz$suicide
w <- df_yz$w_rmpw
e1 <- df_yz$e1
e0 <- df_yz$e0

# OUTCOMES
mu1 <- mean(Y[G == 1])
mu0 <- mean(Y[G == 0])

mu10 <- weighted.mean(Y[G == 1], w[G == 1])

mu1
mu10
mu0

XA <- model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))) %>% NAImpute()
XA <- XA[, !grepl(":.*NA$", colnames(XA))]
XN <- model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% NAImpute()
XN <- XN[, !grepl(":.*NA$", colnames(XN))]


num_cov_lbl = 80; num_pt_lbl = 3; psize = 6; estimand = "reduction"

generatePlot <- function(num_cov_lbl = 8, num_pt_lbl = 3, psize = 6, estimand = "reduction") {
  
  if (estimand == "resid" | estimand == "residual") {
    estimand <- "residual"
    seq <- seq(1, 3, by = 0.01)
    title <- paste0("Disp. RESIDUAL for outcome ", toupper(outcome))
  } else {
    estimand <- "reduction"
    seq <- seq(1, 2, by = 0.01)
    title <- paste0("Disp. REDUCTION for outcome ", toupper(outcome))
  }
  print(title)
  
  out <- sapply(seq, function(i) {
    decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(i), estimand = estimand, RD = TRUE, verbose = FALSE)
  }) %>% t() %>% data.frame() %>% mutate(lam = seq) %>% relocate(lam, .before = "X1")
  
  p_extrema <- out %>% 
    pivot_longer(cols = -lam, names_to = "extrema", values_to = "value") %>%
    ggplot(aes(x = lam, y = value, color = log(lam))) +
    geom_point(size = 1) + 
    theme_minimal()
  
  # find out$X1 where the lower bound crosses 0
  
  # Assuming 'sequence' is your vector of numbers
  
  cross_ind <- max(max(which(out$X1 >= 0)), # last positive index
                   min(which(out$X1 <= 0))) # first negative index
  
  Lam_pe <- out$lam[cross_ind]
  Lam_05 <- 1.02
  print(paste0("Point Est Lambda: ", Lam_pe))
  print(paste0("Point Est Lambda: ", Lam_05))
  # even if we care about red or res, we use point bc the lambda already takes
  # into account whether we care about red or res

  amplification <- decompsens::decompAmplify(G, Z, XA, XN, Y, w, mu_10=mu10, 
                                             Lambda = Lam_pe, e1 = e1, e0 = e0)
  # amp for point est
  bounds <- amplification$maxbias
  
  strongest_cov_df <- amplification[[1]] %>% drop_na()
  max_imbal_stnd <- amplification$max_imbal_stnd
  max_betau_01 <- amplification$max_beta
  
  maxbias <- max(abs(bounds)) 
  max_betau_01 <- max(strongest_cov_df$coeff)
  max_imbal <- max(strongest_cov_df$imbal)
  
  
  # Amp for 0.05
  amplification05 <- decompsens::decompAmplify(G, Z, XA, XN, Y, w, mu_10=mu10, 
                                             Lambda = Lam_05, e1 = e1, e0 = e0)
  maxbias05 <- amplification05$maxbias
  
  if (estimand == "residual") {
    beta <- seq(max(0, min(strongest_cov_df$coeff) - 0.2), max(strongest_cov_df$coeff) + 0.5, by = 0.01)
    imbalance <- seq(0.01, max(strongest_cov_df$imbal, strongest_cov_df$imbal_wt) + 0.2, by = 0.005)
  } else {
    beta <- seq(max(0, min(strongest_cov_df$coeff) - 0.2), max(strongest_cov_df$coeff) + 0.2, by = 0.01)
    imbalance <- seq(0.001, max(strongest_cov_df$imbal, strongest_cov_df$imbal_wt) + 0.2, by = 0.005)
  }
  
  data.fit <- expand.grid(beta, imbalance)
  names(data.fit) <- c("beta", "imbalance")
  
  df_plot <- data.fit %>% 
    mutate(bias = abs(beta * imbalance))
  
  round_down <- function(a) {
    floor(a * 100) / 100
  }
  
  if (estimand == "residual") {
    bins <- seq(round(maxbias, 2) - 0.2, round(maxbias, 2) + 0.1, length.out = 4) %>% round(3)
    bins <- c(0.05, 0.15, 0.2, 0.25)
  } else {
    bins <- seq(round(maxbias, 2) - 0.05, round(maxbias, 2) + 0.05, length.out = 4) %>% round(3)
    bins <- c(0.01, 0.05, 0.1, 0.2)
  }
  
  bins <- bins[bins > 0]
  
  p1 <- df_plot %>%
    filter(beta >= 0) %>% 
    ggplot(aes(x = imbalance, y = beta, z = bias)) +
    geom_contour(col="gray22", breaks = bins, linewidth = 0.6) + 
    metR::geom_text_contour(aes(z = bias), 
                            breaks = bins,
                            stroke = 0.2, skip = 0) + 
    # metR::geom_contour_fill(breaks = c(maxbias05, 1000 * maxbias05), fill='gray55', alpha = 0.2) +
    # geom_contour(breaks = c(maxbias05), col='gray62', linewidth = 1.2, alpha = 2) + 
    # metR::geom_text_contour(aes(z = maxbias05), stroke = 0.2) + 
    metR::geom_contour_fill(breaks = c(maxbias, 1000 * maxbias), fill='dodgerblue3', alpha = switch(estimand, "reduction" = 0.2, "residual" = 0.26)) +
    geom_contour(breaks = c(maxbias), col='dodgerblue3', linewidth = 1.2, alpha = 2) + 
    metR::geom_text_contour(aes(z = maxbias), stroke = 0.2) + theme_minimal()
  #geom_point(x = 0.65, y = maxbias / 0.65, size = psize - 2, color = "black")
  
  # geom_contour(data = data.frame(x = seq(0.01, 0.25, by = 0.01), y = seq(0.01, 2, by = 0.01)) %>% 
  #                mutate(z = x * y),
  #              aes(x = x, y = y, z = z), breaks = c(maxbias), col='blue', linewidth = 1, inherit.aes = F) +
  
  p1
  
  strongest_cov_df_long <- strongest_cov_df %>% 
    pivot_longer(cols = c("imbal", "imbal_wt"), names_to = "imbal_type", values_to = "imbal_val") #%>% filter(imbal_type == "imbal_wt")
  
  num_cov <- min(nrow(strongest_cov_df), num_cov_lbl) * 2
  p1_full <- p1 + geom_point(data = (strongest_cov_df_long %>% filter(covar %in% non_allowable_covs) %>% slice(1:num_cov)), 
                             aes(x = imbal_val, y = coeff, z = 0, color = imbal_type), size = psize) + 
    scale_color_manual(labels = c("imbal" = "Pre-wt", "imbal_wt" = "Post-wt"),
                       values = c("imbal" = "red", "imbal_wt" = "forestgreen"))
  p1_full
  
  # teal: #00BFC4
  # reddish: #F8766D
  
  hull <- strongest_cov_df %>%
    slice(chull(coeff, imbal))
  
  p1_full_unscaled <- p1_full + 
    # ggrepel::geom_label_repel(data = strongest_cov_df %>% filter(covar %in% non_allowable_covs) %>% slice(1:num_pt_lbl), 
    #                           aes(x = imbal_wt, y = coeff, z = 0, label = covar),
    #                           nudge_y = 0.015, nudge_x = 0.005, label.padding = 0.05,
    #                           point.padding = 0.05, size = 5.5) + 
    ggrepel::geom_label_repel(data = strongest_cov_df %>% filter(covar %in% non_allowable_covs) %>% slice(1:num_pt_lbl), 
                              aes(x = imbal, y = coeff, z = 0, label = covar),
                              nudge_y = switch(estimand, "reduction" = 0.015, "residual" = 0.03), 
                              nudge_x = switch(estimand, "reduction" = 0.005, "residual" = 0.01), 
                              label.padding = 0.05,
                              point.padding = 0.05, size = 5.5) + 
    # ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov_lbl,],
    #                           aes(x = imbal_wt, y = coeff, z = 0, label = covar),
    #                           nudge_y = 0.005, nudge_x = 0.005, label.padding = 0.1,
    #                           point.padding = 0.1) + 
    #geom_polygon(data = hull, aes(x = imbal, y = coeff, z = 0), alpha = 0.3, fill = "red") + 
    theme_minimal(base_size = 30) + 
    #theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
    labs(x = TeX("$|\\delta_u|$ Imbalance"), y = TeX("$|\\beta_u|$ Impact"), 
         color = c("Imbalance")) + 
    #ylim(0, 1) + 
    geom_hline(yintercept = 0, color = "gray55") + scale_y_continuous(limits = c(0, NA)) + 
    theme(legend.position="bottom")
  #+ xlim(0, 1) + ylim(0, 1)
  #+ scale_y_continuous(limits = c(0, NA))
  
  if (estimand == "reduction") {
    p_out <- p1_full_unscaled + 
      geom_text(x = 0.65, y = 0.23,
                label = TeX(paste0("$\\Lambda^{*} = ", Lam_pe)),
                size = 8, color = "black") +
      geom_text(x = 0.65, y = 0.2,
                label = paste0("Bias: ", round(maxbias, 3)), size = 8, color = "black") 
    p_out
  } else {
    p_out <- p1_full_unscaled + 
      geom_text(x = 0.7, y = 0.53,
                label = TeX(paste0("$\\Lambda^{*} = ", Lam_pe)),
                size = 8, color = "black") +
      geom_text(x = 0.7, y = 0.48,
                label = paste0("Bias: ", round(maxbias, 3)), size = 8, color = "black") 
    p_out
  }
}

red_plot <- generatePlot(num_cov_lbl = 12, num_pt_lbl = 3, psize = 6.5, estimand = "red")
red_plot

resid_plot <- generatePlot(num_cov_lbl = 12, num_pt_lbl = 3, psize = 6.5, estimand = "resid")
resid_plot

