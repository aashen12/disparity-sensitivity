# Analysis of list 1
library(tidyverse)
library(decompsens)
library(parallel)
library(doParallel)
library(boot)
library(latex2exp)

source("functions.R")
set.seed(122357)

numCores <- parallel::detectCores()

doParallel::registerDoParallel(numCores)

Z_method <- "better_upset"
# "better_worry",
# "smile",
# "better_upset",
# "love",
# "easy_talk
# aggregate



df_x <- read_csv(paste0("../data/list1_X_", Z_method, ".csv"))
df_yz <- read_csv(paste0("../data/list1_YZG_", Z_method, ".csv"))

G <- df_yz$sex_min
Z <- df_yz$parent_accept
Y <- df_yz$ideation

# Allowable covariates 
options(na.action='na.pass')


#mediators <- c("peer_victimization")
mediators <- c("")
allowable_covs <- c("age", "sex", "sib_num", "sib_order")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df <- df_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute()
XA_log <- create_model_matrix(df_allowable)

df <- df_non_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute()
XN_log <- create_model_matrix(df_non_allowable)

# log stands for "logistic". We use a more comlex martix to construct the propensity scores

# XA_log <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>%
#   data.frame() %>%
#   NAImpute()
# XA_log <- XA_log[, !grepl(":.*NA$", colnames(XA_log))]
# XN_log <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% NAImpute()
# XN_log <- XN_log[, !grepl(":.*NA$", colnames(XN_log))]

X_log <- cbind(XA_log, XN_log)

# Non-allowable covariates


df_yz %>% group_by(sex_min) %>% 
  summarise(
    n = n(),
    mean_parent_accept = mean(parent_accept, na.rm = FALSE),
    sd_parent_accept = sd(parent_accept, na.rm = FALSE)
  )

allowable <- TRUE
trim <- 0.01

if (allowable) {
  if (Z_method == "aggregate") {
    e0 <- glm(Z ~ ., data = XA_log, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ ., data = X_log, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  } else {
    e0 <- glm(Z ~ XA_log, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ X_log, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  }
} else {
  e0 <- glm(Z ~ X_log, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
  e1 <- glm(Z ~ X_log, family = binomial, weights = G, na.action = na.exclude)$fitted.values
}

e0 <- pmax(pmin(e0, 1 - trim), trim)
e1 <- pmax(pmin(e1, 1 - trim), trim)

#wr <- glm(G ~ 1, family = binomial, na.action = na.exclude)$fitted.values / glm(G ~ XA[, "age"] + XA[, "sexF"] + XA[, "sexM"], family = binomial, na.action = na.exclude)$fitted.values
wr <- 1

w1 = e0 / e1 * wr
w0 = (1 - e0) / (1 - e1) * wr
w = w1 * Z + w0 * (1 - Z)

summary(w)

df_yz <- df_yz %>% mutate(w_rmpw = w)
write_csv(df_yz, paste0("../data/list1_YZGW_", Z_method, ".csv"))

message(paste0("CSV file for Z method ", Z_method, " has been written."))

# assess balance in X



X_plot <- X_plot <- cbind(model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))),
                          model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs)))) %>% NAImpute()
X_stnd <- apply(X_plot, 2, scale)

pre_weight_G <- colMeans(X_stnd[G == 1, ]) - colMeans(X_stnd[G == 0, ])
max(abs(pre_weight_G))

post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - 
  colSums(X_stnd[G == 0, ] * w[G == 0])  / sum(w[G == 0])
lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")

post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - colMeans(X_stnd[G == 1, ])
lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")


pre_weight_Z <- colMeans(X_stnd[Z == 1, ]) - colMeans(X_stnd[Z == 0, ])
max(abs(pre_weight_Z))

post_weight_Z <- colSums(X_stnd[Z == 1, ] * w[Z == 1] / sum(w[Z == 1])) - colMeans(X_stnd[Z == 0, ])
max(abs(post_weight_Z))

loveZ <- lovePlot(pre_weight_Z, post_weight_Z, title = "Covariate Balance wrt Parental Support")
loveZ

# Assess balance in X for the sexual minority group

# pre_weight_G1 <- colMeans(X_stnd[G == 1 & Z == 1, ]) - colMeans(X_stnd[G == 1 & Z == 0, ])
# max(abs(pre_weight_G1))
# 
# post_weight_G1 <- colSums(X_stnd[G == 1 & Z == 1, ] * w[G == 1 & Z == 1] / sum(w[G == 1 & Z == 1])) - 
#   colSums(X_stnd[G == 1 & Z == 0, ] * w[G == 1 & Z == 0] / sum(w[G == 1 & Z == 0]))
# max(abs(post_weight_G1))
# 
# loveG1 <- lovePlot(pre_weight_G1, post_weight_G1, title = "Covariate Balance wrt Parental Support for Sex. Minorities ONLY")
# loveG1
# 
# pre_weight_G0 <- colMeans(X_stnd[G == 0 & Z == 1, ]) - colMeans(X_stnd[G == 0 & Z == 0, ])
# max(abs(pre_weight_G0))
# 
# post_weight_G0 <- colSums(X_stnd[G == 0 & Z == 1, ] * w[G == 0 & Z == 1] / sum(w[G == 0 & Z == 1])) - 
#   colSums(X_stnd[G == 0 & Z == 0, ] * w[G == 0 & Z == 0] / sum(w[G == 0 & Z == 0]))
# max(abs(post_weight_G0))
# 
# loveG0 <- lovePlot(pre_weight_G0, post_weight_G0, title = "Covariate Balance wrt Parental Support for non sex. minorities")
# loveG0
# 
# pre_weight_Z1 <- colMeans(X_stnd[Z == 1 & G == 1, ]) - colMeans(X_stnd[Z == 1 & G == 0, ])
# max(abs(pre_weight_Z1))
# 
# post_weight_Z1 <- colSums(X_stnd[Z == 1 & G == 1, ] * w[Z == 1 & G == 1] / sum(w[Z == 1 & G == 1])) - 
#   colSums(X_stnd[Z == 1 & G == 0, ] * w[Z == 1 & G == 0] / sum(w[Z == 1 & G == 0]))
# max(abs(post_weight_Z1))
# 
# loveZ1 <- lovePlot(pre_weight_Z1, post_weight_Z1, title = "Covariate Balance wrt sex min status for treated individuals")
# loveZ1
# 
# pre_weight_Z0 <- colMeans(X_stnd[Z == 0 & G == 1, ]) - colMeans(X_stnd[Z == 0 & G == 0, ])
# max(abs(pre_weight_Z0))
# 
# post_weight_Z0 <- colSums(X_stnd[Z == 0 & G == 1, ] * w[Z == 0 & G == 1] / sum(w[Z == 0 & G == 1])) - 
#   colSums(X_stnd[Z == 0 & G == 0, ] * w[Z == 0 & G == 0] / sum(w[Z == 0 & G == 0]))
# max(abs(post_weight_Z0))
# 
# loveZ0 <- lovePlot(pre_weight_Z0, post_weight_Z0, title = "Covariate Balance wrt sex min status for control individuals")
# loveZ0
# 
# 
# pre_weight_subopt <- colMeans(X_stnd[Z == 0 & G == 1, ]) - colMeans(X_stnd[Z == 1 & G == 0, ])
# 
# post_weight_subopt <- colSums(X_stnd[Z == 0 & G == 1, ] * w[Z == 0 & G == 1] / sum(w[Z == 0 & G == 1])) - 
#   colSums(X_stnd[Z == 1 & G == 0, ] * w[Z == 1 & G == 0] / sum(w[Z == 1 & G == 0]))
# 
# lovesubopt <- lovePlot(pre_weight_subopt, post_weight_subopt, title = "Covariate Balance for the suboptimal situation")
# lovesubopt





