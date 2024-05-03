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
allowable_covs <- c("age", "sex", "sib_num", "sib_order")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df <- df_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

df <- df_non_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

interact <- FALSE

if (interact == TRUE) {
  XA_log <- create_model_matrix(df_allowable)
  XN_log <- create_model_matrix(df_non_allowable)
  X_log <- cbind(XA_log, XN_log)
} else {
  XA_log <- df_allowable
  XN_log <- df_non_allowable
  X_log <- cbind(XA_log, XN_log)
}


G <- df_yz$sex_min
Z <- df_yz$parent_accept
Y <- df_yz$suicide
w <- df_yz$w_rmpw
e1 <- df_yz$e1
e0 <- df_yz$e0

# OUTCOMES
mu1 <- mean(Y[G == 1])
mu0 <- mean(Y[G == 0])

mean(Z[G == 1]) - mean(Z[G == 0])

mean(Z[G == 1]) / mean(Z[G == 0])




mu1
mu0
obs_disp <- mu1 - mu0
obs_disp
mean(Y[G == 1]) - mean(Y[G == 0])

sd_obs_disp <- sqrt(var(Y[G == 1]) / sum(G == 1) + var(Y[G == 0]) / sum(G == 0))
# 
# construct 95% CI for obs_disp
ci <- c(obs_disp - qnorm(0.975) * sd_obs_disp, obs_disp + qnorm(0.975) * sd_obs_disp)
ci


mu10 <- weighted.mean(Y[G == 1], w[G == 1])


mu1_AS <- lm(Y ~ df_allowable$age + df_allowable$sexF, weights = G)$fitted.values %>% mean()
mu0_AS <- lm(Y ~ df_allowable$age + df_allowable$sexF, weights = 1-G)$fitted.values %>% mean()

mu1_AS
mu0_AS


mu1
mu10
mu0



obs_disp

reduction <- mu1 - mu10
reduction
mu1/mu10
reduction / obs_disp


residual <- mu10 - mu0
residual
mu10/mu0
residual / obs_disp


# Bootstrap standard errors
# B <- 1000
# resid_boot <- numeric(B)
# allowable <- TRUE
# out <- parallel::mclapply(1:B, function(i) {
#   ind <- sample(1:length(Y), length(Y), replace = TRUE)
#   w_boot <- decompsens::estimateRMPW(G=G[ind], Z=Z[ind], Y=Y[ind], XA=XA_log[ind,], XN=XN_log[ind,],
#                                      trim = 0.01, allowable = TRUE)
# 
#   mu10_boot <- sum(Y[ind][G[ind] == 1] * w_boot[G[ind] == 1]) / sum(w_boot[G[ind] == 1])
#   mu1_boot <- mean(Y[ind][G[ind] == 1])
#   mu0_boot <- mean(Y[ind][G[ind] == 0])
#   resid_boot <- mu10_boot - mu0_boot
#   red_boot <- mu1_boot - mu10_boot
#   list(resid_boot = resid_boot, red_boot = red_boot)
# }, mc.cores = numCores)
# 
# out_resid <- unlist(lapply(out, function(x) x[["resid_boot"]]))
# out_red <- unlist(lapply(out, function(x) x[["red_boot"]]))
# 
# quantile(out_resid, c(0.025, 0.975))
# c(2*residual - quantile(out_resid, c(0.975, 0.025)), residual) # pivot bootstrap CI
# 
# quantile(out_red, c(0.025, 0.975))
# c(2*reduction - quantile(out_red, c(0.975, 0.025)), reduction)
# 
# sd(out_red)
# mean(out_red)



seq <- seq(1, 3, by = 0.005)
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
  ggplot(aes(x = lam, y = value, color = lam)) +
  geom_point(size = 2) + 
  facet_wrap(~extrema) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = mu0, color = "green", linewidth = 1) +
  geom_hline(yintercept = mu1, color = "green", linewidth = 1) +
  geom_hline(yintercept = mu10, color = "red", linewidth = 1) +
  theme_minimal(base_size = 25)


# make one for point estimate!!

p_extrema


decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(1.12), estimand = "red", RD = TRUE, verbose = FALSE)
mu1 - mu10

decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(1.05), estimand = "point", RD = TRUE, verbose = FALSE)
mu1
