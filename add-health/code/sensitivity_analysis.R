# perform analysis of add health data
# devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))

library(tidyverse)
library(decompsens)
library(latex2exp)

set.seed(122)

Hajek <- function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}

HT <- function(x, w, na.rm = TRUE) {
  mean(x * w, na.rm = na.rm)
}


df_yz <- read_csv("../../data/add-health-data/trt_out_df.csv") %>% select(-attempt) %>% drop_na()
df_x <- read_csv("../../data/add-health-data/covariate_df.csv") %>% select(-dad_ed) %>% drop_na()

covariates <- names(df_x)[!names(df_x) %in% c("AID", "CLUSTER2", "GSWGT1")]

allowables <- c("age", "sex")

non_allowables <- setdiff(covariates, allowables)

df_full <- inner_join(df_yz, df_x, by = "AID")
# G_tot <- df_full$sex_minority
# G1_ind <- which(G_tot == 1)
# G0_ind <- sample(which(G_tot == 0), size = sum(G_tot == 1), replace = TRUE)
# 
# df <- df_full %>% slice(c(G1_ind, G0_ind))

df <- df_full

XA <- model.matrix(~ .^2 -1, data = df %>% select(all_of(allowables)))
XN <- model.matrix(~ .^2 -1, data = df %>% select(all_of(non_allowables)))

X <- cbind(XA, XN)

head(X)

Y <- df$ideation
Z <- df$school_belong
G <- df$sex_minority

mu1 <- mean(Y[G == 1])
mu0 <- mean(Y[G == 0])



obs_disp <- mu1 - mu0
obs_disp

sd_obs_disp <- sqrt(var(Y[G == 1]) / sum(G == 1) + var(Y[G == 0]) / sum(G == 0))

# construct 95% CI for obs_disp
alpha <- 0.05
ci <- c(obs_disp - qnorm(1-alpha/2) * sd_obs_disp, obs_disp + qnorm(1-alpha/2) * sd_obs_disp)
ci


mu1/mu0

trim <- 0
allow <- T

w <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=XA, XN=XN, trim = trim, allowable = allow)
summary(w)

mu10 <- sum(Y[G == 1] * w[G == 1]) / sum(w[G == 1])
mu10

mu1
mu0

mu1-mu0


reduction <- mu1 - mu10
reduction
mu1/mu10

residual <- mu10 - mu0
residual
mu10/mu0

# ## Bootstrap standard errors
# B <- 2000
# resid_boot <- numeric(B)
# out <- parallel::mclapply(1:B, function(i) {
#   ind <- sample(1:length(Y), length(Y), replace = TRUE)
#   w_boot <- decompsens::estimateRMPW(G=G[ind], Z=Z[ind], Y=Y[ind], XA=XA[ind,], XN=XN[ind,], trim = trim, allowable = allow)
#   mu10_boot <- sum(Y[ind][G[ind] == 1] * w_boot[G[ind] == 1]) / sum(w_boot[G[ind] == 1])
#   mu1_boot <- mean(Y[ind][G[ind] == 1])
#   mu0_boot <- mean(Y[ind][G[ind] == 0])
#   resid_boot <- mu10_boot - mu0_boot
#   red_boot <- mu1_boot - mu10_boot
#   list(resid_boot = resid_boot, red_boot = red_boot)
# })
# 
# out_resid <- unlist(lapply(out, function(x) x[["resid_boot"]]))
# out_red <- unlist(lapply(out, function(x) x[["red_boot"]]))
# 
# se_resid <- sqrt((1 / (B-1)) * sum((out_resid - mean(out_resid))^2))
# se_resid
# # 90% CI
# ci_resid <- c(residual - qnorm(0.95) * se_resid, residual + qnorm(0.95) * se_resid)
# ci_resid
# 
# se_red <- sqrt((1 / (B-1)) * sum((out_red - mean(out_red))^2))
# se_red
# # 90% CI
# ci_red <- c(reduction - qnorm(0.8) * se_red, reduction + qnorm(0.8) * se_red)
# ci_red


###########################################################################################
###########################################################################################
###########################################################################################
lam <- 5
h <- log(lam)

wh <- w * exp(h)

mu10_h <- weighted.mean(Y[G == 1], wh[G == 1])
mu10_h
mu10

mu0
mu1


seq <- seq(1, 2, by = 0.005)
out <- sapply(seq, function(i) {
  red <- decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(i), estimand = "red", RD = TRUE, verbose = FALSE)
  res <- decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(i), estimand = "res", RD = TRUE, verbose = FALSE)
  point <- decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(i), estimand = "point", RD = TRUE, verbose = FALSE)
  c(red = red, res = res, point = point)
}) %>% t() %>% data.frame() %>% mutate(lam = seq, .before = "red1")

p_extrema <- out %>% 
  pivot_longer(cols = -lam, names_to = "extrema", values_to = "value") %>%
  mutate(extrema = gsub("[0-9]", "", extrema)) %>%
  mutate(point_est = case_when(
    extrema == "red" ~ mu1 - mu10,
    extrema == "res" ~ mu10 - mu0,
    extrema == "point" ~ mu10
  )) %>% 
  ggplot(aes(x = lam, y = value, color = log(lam))) +
  geom_point(size = 2) + 
  facet_wrap(~extrema) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = mu0, color = "green", linewidth = 1) +
  geom_hline(yintercept = mu1, color = "green", linewidth = 1) +
  geom_hline(yintercept = mu10, color = "red", linewidth = 1) +
  theme_minimal(base_size = 25)


# make one for point estimate!!

p_extrema


decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(1.06), estimand = "red", RD = TRUE, verbose = FALSE)
mu1 - mu10

decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(1.05), estimand = "point", RD = TRUE, verbose = FALSE)
mu1

num_cov_lbl <- 4

psize <- 6


generatePlot <- function(num_cov_lbl = 8, psize = 6, estimand = "resid") {
  
  if (estimand == "resid") {
    estimand <- "residual"
    seq <- seq(1, 3, by = 0.01)
    title <- "Disp. RESIDUAL (ADD)"
  } else {
    estimand <- "reduction"
    seq <- seq(1, 2, by = 0.01)
    title <- "Disp. REDUCTION (ADD)"
  }
  
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
  
  cross_ind <- min(max(which(out$X1 >= 0)), # last positive index
                   min(which(out$X1 <= 0))) # first negative index
  
  Lam <- out$lam[cross_ind]
  
  print(paste0("Crit Lam for ", estimand, ": ",Lam))
  
  bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, Lambda = Lam)
  amplification <- decompsens::informalAmplify(G, Z, XA, XN, Y, Lambda = Lam)
  
  strongest_cov_df <- amplification[[1]] %>% drop_na()
  max_imbal_stnd <- amplification$max_imbal_stnd
  max_betau_01 <- amplification$max_beta
  
  maxbias <- max(abs(bounds[[1]]))
  
  max_betau_01 <- max(strongest_cov_df$coeff)
  max_imbal <- max(strongest_cov_df$imbal)
  
  if (estimand == "residual") {
    beta <- seq(min(strongest_cov_df$coeff) - 0.2, max(strongest_cov_df$coeff) + 0.2, by = 0.01)
    imbalance <- seq(0.001, max(strongest_cov_df$imbal) + 0.2, by = 0.005)
  } else {
    beta <- seq(min(strongest_cov_df$coeff) - 0.2, max(strongest_cov_df$coeff) + 0.2, by = 0.01)
    imbalance <- seq(0.001, max(strongest_cov_df$imbal) + 0.2, by = 0.005)
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
    bins <- c(0.05, 0.15, 0.3)
  } else {
    bins <- seq(round(maxbias, 2) - 0.05, round(maxbias, 2) + 0.05, length.out = 4) %>% round(3)
    bins <- c(0.05, 0.12, 0.22)
  }
  
  bins <- bins[bins > 0]
  
  
  
  p1 <- df_plot %>%
    ggplot(aes(x = imbalance, y = beta, z = bias)) +
    theme_bw() + 
    geom_contour(col="gray55", breaks = bins) + 
    metR::geom_text_contour(aes(z = bias), 
                            breaks = bins,
                            stroke = 0.2, skip = 0) + 
    metR::geom_contour_fill(breaks = c(maxbias, 1000 * maxbias), fill='powderblue', alpha = 0.5) +
    geom_contour(breaks = c(maxbias), col='blue', linewidth = 1) + 
    metR::geom_text_contour(aes(z = maxbias), stroke = 0.2) + 
    geom_point(x = 0.65, y = maxbias / 0.65, size = psize - 2, color = "black")
    
  
  # geom_contour(data = data.frame(x = seq(0.01, 0.25, by = 0.01), y = seq(0.01, 2, by = 0.01)) %>% 
  #                mutate(z = x * y),
  #              aes(x = x, y = y, z = z), breaks = c(maxbias), col='blue', linewidth = 1, inherit.aes = F) +
  
  p1
  
  strongest_cov_df_long <- strongest_cov_df %>% 
    pivot_longer(cols = c("imbal", "imbal_wt"), names_to = "imbal_type", values_to = "imbal_val") 
  
  num_cov <- max(nrow(strongest_cov_df), num_cov_lbl)
  p1_full <- p1 + geom_point(data = strongest_cov_df_long %>% filter(imbal_type == "imbal") %>% slice(1:num_cov), 
                             aes(x = imbal_val, y = coeff, z = 0), size = psize, color = "red")
  
  # teal: #00BFC4
  # reddish: #F8766D
  
  hull <- strongest_cov_df[1:(num_cov/2),] %>%
    slice(chull(coeff, imbal))
  
  
  
  p1_full_unscaled <- p1_full + 
    ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov_lbl,], 
                              aes(x = imbal, y = coeff, z = 0, label = covar),
                              nudge_y = -0.005, nudge_x = 0.005, label.padding = 0.1,
                              point.padding = 0.1) + 
    # ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov_lbl,],
    #                           aes(x = imbal_wt, y = coeff, z = 0, label = covar),
    #                           nudge_y = 0.005, nudge_x = 0.005, label.padding = 0.1,
    #                           point.padding = 0.1) + 
    geom_polygon(data = hull, aes(x = imbal, y = coeff, z = 0), alpha = 0.3, fill = "red") + 
    theme_bw(base_size = 20) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
    labs(x = TeX("absolute $\\delta_u$"), y = TeX("$\\beta_u$"),
         title = paste0(title), color = c("Imbalance")) + 
    theme(legend.position = "bottom") + 
    geom_hline(yintercept = 0, color = "gray55") + 
    geom_text(x = 0.55 - 0.1, y = 0.365,
              label = TeX(paste0("$\\Lambda^{*} = ", Lam)), 
              size = 8, color = "black") + 
    geom_text(x = 0.55 - 0.1, y = 0.365 - 0.04,
              label = paste0("Bias: ", round(maxbias, 3)), size = 8, color = "black")
  
  p1_full_unscaled
}

resid_plot <- generatePlot(num_cov_lbl = 8, psize = 5, estimand = "resid")
red_plot <- generatePlot(num_cov_lbl = 8, psize = 5, estimand = "red")

resid_plot
ggsave("plots/resid_plot_add.png", resid_plot, width = 10, height = 10)

red_plot
ggsave("plots/red_plot_add.png", red_plot, width = 10, height = 10)

