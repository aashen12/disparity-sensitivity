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
#allowable_covs <- c("age", "sex", "sib_order", "sib_num")
allowable_covs <- c("age", "sex")
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
B <- 1000
resid_boot <- numeric(B)
allowable <- TRUE

mu10
getExtrema(G,Y,gamma = log(1), w = w, verbose = FALSE, estimand = "point")


reduction
getExtrema(G,Y,gamma = log(1.09), w = w, verbose = TRUE, estimand = "red")

residual
getExtrema(G,Y,gamma = log(1.67), w = w, verbose = TRUE, estimand = "res")


e0_raw <- glm(Z ~ ., data = df_allowable, family = binomial(link = "logit"), weights = (1-G))$fitted.values
e1_raw <- glm(Z ~ ., data = cbind(df_allowable, df_non_allowable), family = binomial(link = "logit"), weights = G)$fitted.values


trim0 <- 1 - quantile(e0_raw, 0.99, names = FALSE)
trim1 <- switch(outcome, 
                "ideation" = 1 - quantile(e1_raw, 0.93, names = FALSE),
                "attempt" = 1 - quantile(e1_raw, 0.90, names = FALSE))


boot_ci_red <- bootstrapCI(G, Z, Y, XA_log, XN_log,
                           gamma = log(1.02), trim0 = trim0, trim1 = trim1,
                           estimand = "red", stratify = FALSE, alpha = 0.05,
                           allowable = TRUE)
boot_ci_red

# Want to find Lambda such that mu10 = mu0. # 1.3 for ideation, 1.7 for attempt
mu10
mu0
bootstrapCI(G, Z, Y, XA_log, XN_log, gamma = log(1), 
            trim0 = trim0, trim1 = trim1,
            estimand = "point", alpha = 0.1)




out <- parallel::mclapply(1:B, function(i) {
  ind <- sample(1:length(Y), length(Y), replace = TRUE)
  # ind_G1 <- sample(which(G == 1), sum(G == 1), replace = TRUE)
  # ind_G0 <- sample(which(G == 0), sum(G == 0), replace = TRUE)
  # ind <- c(ind_G1, ind_G0)
  w_boot_obj <- decompsens::estimateRMPW(G=G[ind], Z=Z[ind], Y=Y[ind],
                                         XA=XA_log[ind,], XN=XN_log[ind,],
                                         trim0 = trim0, trim1 = trim1,
                                         allowable = TRUE)
  w_boot <- w_boot_obj$w_rmpw
  # w_boot <- w[ind]
  mu10_boot <- sum(Y[ind][G[ind] == 1] * w_boot[G[ind] == 1]) / sum(w_boot[G[ind] == 1])
  mu1_boot <- mean(Y[ind][G[ind] == 1])
  mu0_boot <- mean(Y[ind][G[ind] == 0])
  resid_boot <- mu10_boot - mu0_boot
  red_boot <- mu1_boot - mu10_boot
  list(resid_boot = resid_boot, red_boot = red_boot)
}, mc.cores = numCores)

out_resid <- unlist(lapply(out, function(x) x[["resid_boot"]]))
out_red <- unlist(lapply(out, function(x) x[["red_boot"]]))

quantile(out_red, c(0.025, 0.975)) # percentile bootstrap
reduction
mean(out_red) + c(-1, 1) * qnorm(0.975) * sd(out_red) # boot CI
c(2*reduction - quantile(out_red, c(0.975, 0.025)), reduction) # pivot bootstrap

quantile(out_resid, c(0.025, 0.975)) # percentile bootstrap
mean(out_resid) + c(-1, 1) * qnorm(0.975) * sd(out_resid) # boot CI
c(2*residual - quantile(out_resid, c(0.975, 0.025)), residual) # pivot bootstrap


reduction
sd(out_red)
mean(out_red)

# 95% CI for out_red



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

# 05/19
# make one for point estimate!!

p_extrema

