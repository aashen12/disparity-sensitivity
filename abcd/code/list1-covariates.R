# Analysis of list 1
rm(list = ls())

library(tidyverse)
library(decompsens)
library(parallel)
library(doParallel)
library(boot)
library(cobalt)
library(latex2exp)
library(WeightIt)
library(lsr)

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

outcome <- "attempt" # ideation or attempt


interact <- FALSE

df_x <- read_csv(paste0("../data/list1_X_", Z_method, "_", outcome, ".csv"))
df_yz <- read_csv(paste0("../data/list1_YZG_", Z_method, "_", outcome, ".csv"))

G <- df_yz$sex_min
Z <- df_yz$parent_accept
Y <- df_yz$suicide

table(Z)
table(G)
table(Y)

# Allowable covariates 
options(na.action='na.pass')


#mediators <- c("peer_victimization")
mediators <- c("src_subject_id", "family_mental_health")
allowable_covs <- c("age", "sex", "sib_order", "sib_num")
# allowable_covs <- c("age", "sex")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

df_non_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

if (interact == TRUE) {
  XA_log <- create_model_matrix(df_allowable)
  XN_log <- create_model_matrix(df_non_allowable)
  X_log <- cbind(XA_log, XN_log)
}

# log stands for "logistic". We use a more comlex martix to construct the propensity scores

# XA_log <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>%
#   data.frame() %>%
#   NAImpute()
# XA_log <- XA_log[, !grepl(":.*NA$", colnames(XA_log))]
# XN_log <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% NAImpute()
# XN_log <- XN_log[, !grepl(":.*NA$", colnames(XN_log))]


# Non-allowable covariates


df_yz %>% group_by(sex_min) %>% 
  summarise(
    n = n(),
    mean_parent_accept = mean(parent_accept, na.rm = FALSE),
    sd_parent_accept = sd(parent_accept, na.rm = FALSE)
  )

X_plot <- cbind(model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))),
                model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs)))) %>% NAImpute()

allow <- TRUE

e0_raw <- glm(Z ~ ., data = df_allowable, family = binomial(link = "logit"), weights = (1-G))$fitted.values
e1_raw <- glm(Z ~ ., data = cbind(df_allowable, df_non_allowable), family = binomial(link = "logit"), weights = G)$fitted.values


trim0 <- 1 - quantile(e0_raw, 0.99, names = FALSE)
trim1 <- switch(outcome, 
                "ideation" = 1 - quantile(e1_raw, 0.93, names = FALSE),
                "attempt" = 1 - quantile(e1_raw, 0.90, names = FALSE))


e0 <- pmax(pmin(e0_raw, 1 - trim0), trim0)
e1 <- pmax(pmin(e1_raw, 1 - trim1), trim1)


weighted.var <- function(x, w) {
  #w <- abs(w)
  xbar_w <- weighted.mean(x, w)
  coef <- sum(w) / ( (sum(w))^2 - sum(w^2) )
  coef * sum(w * (x - xbar_w)^2)
}

## For e_0 ###

pre_weight_e0 <- apply(df_allowable, 2, function(col) {
  abs(cohensD(col[G == 0 & Z == 1], col[G == 0 & Z == 0]))
})

post_weight_e0 <- apply(df_allowable, 2, function(col) {
  num <- weighted.mean(col[G == 0 & Z == 1], 1/(e0[G == 0 & Z == 1])) - 
    weighted.mean(col[G == 0 & Z == 0], 1/(e0[G == 0 & Z == 0]))
  den <- sqrt(weighted.var(col[G == 0 & Z == 1], 1/(e0[G == 0 & Z == 1])) + 
                weighted.var(col[G == 0 & Z == 0], 1/(e0[G == 0 & Z == 0]))) * sqrt(0.5)
  abs(num) / den
})

lovePlot(pre_weight_e0, post_weight_e0)


### For e_1 ###

pre_weight_e1 <- apply(X_plot, 2, function(col) {
  abs(cohensD(col[G == 1 & Z == 1], col[G == 1 & Z == 0]))
})

post_weight_e1 <- apply(X_plot, 2, function(col) {
  num <- weighted.mean(col[G == 1 & Z == 1], 1/e1[G == 1 & Z == 1]) - weighted.mean(col[G == 1 & Z == 0], 1/e1[G == 1 & Z == 0])
  den <- sqrt(weighted.var(col[G == 1 & Z == 1], 1/e1[G == 1 & Z == 1]) + 
                weighted.var(col[G == 1 & Z == 0], 1/e1[G == 1 & Z == 0])) / sqrt(2)
  abs(num) / den
})

lovePlot(pre_weight_e1, post_weight_e1)

summary(e1)

w1 = e0 / e1
w0 = (1 - e0) / (1 - e1)
w_rmpw = w1 * Z + w0 * (1 - Z)

summary(w_rmpw[G == 1 & Y == 1])

df_yz <- df_yz %>% mutate(w_rmpw = w_rmpw,
                          e1 = e1,
                          e0 = e0)
write_csv(df_yz, paste0("../data/list1_YZGW_", Z_method, "_", outcome, ".csv"))
message(paste0("CSV file for Z method ", Z_method, " and outcome ", outcome, " has been written."))


# should be equiv
e0_noby <- weightit(Z ~ ., data = df_allowable, method = "glm") #%>% trim(at = 2, lower = TRUE)
e0_test <- glm(Z ~ ., data = df_allowable, family = binomial(link = "logit"))$fitted.values
summary(e0_noby$ps - e0_test)


##################################################################
##################################################################
##################################################################


# Using WeightIt and Cobalt

set.cobalt.options(binary = "std")

df_e0 <- cbind(Z = Z, by = G, df_allowable) %>% filter(by == 0)
df_e1 <- cbind(Z = Z, by = G, df_non_allowable, df_allowable) %>% filter(by == 1)

e0_obj <- weightit(Z ~ . - by, data = df_e0, method = "ebal") #%>% trim(at = 2, lower = TRUE)
# by computes propensity scores within each group of by; different from weights in glm()
e1_obj <- weightit(Z ~ . - by,
                   data = df_e1, method = "ebal") #%>% trim(at = 25, lower = TRUE)






love.plot(e1_obj, drop.distance = TRUE, 
          var.order = "unadjusted",
          abs = TRUE, line = TRUE,
          thresholds = c(m = 0.1), size = 4)


love.plot(e0_obj, drop.distance = TRUE, 
          var.order = "unadjusted",
          abs = TRUE, line = TRUE,
          thresholds = c(m = 0.1),
          size = 6)

summary(e0_obj)
summary(e0_noby)
summary(e1_obj)

bal.tab(e0_obj, un = TRUE, thresholds = c(m = 0.1))
bal.tab(e1_obj, un = TRUE, thresholds = c(m = 0.1))

bal.plot(e1_obj)

e0 <- 1 / e0_obj$weights
e1 <- 1 / e1_obj$weights

summary(e1)

w_rmpw <- (e0 / e1) * Z + ((1-e0) / (1-e1)) * (1-Z)


summary(w_rmpw[G == 1])

# df_yz <- df_yz %>% mutate(w_rmpw = w_rmpw,
#                           e1 = e1,
#                           e0 = e0)
# write_csv(df_yz, paste0("../data/list1_YZGW_", Z_method, "_", outcome, ".csv"))
# message(paste0("CSV file for Z method ", Z_method, " and outcome ", outcome, " has been written."))



# weight_object <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=df_allowable, XN=df_non_allowable,
#                               trim = switch(outcome, "ideation" = 0.05, "attempt" = 0.05), 
#                               allowable = allow)
# 
# 
# 
# w <- weight_object$w_rmpw
# e1 <- weight_object$e1
# e0 <- weight_object$e0

# w1 <- glm(G ~ 1, data = df_allowable, family = binomial(link = "logit"))$fitted.values / 
#   glm(G ~ ., data = df_allowable, family = binomial(link = "logit"))$fitted.values
# w0 <- glm(1 - G ~ 1, data = df_allowable, family = binomial(link = "logit"))$fitted.values / 
#   glm(1 - G ~ ., data = df_allowable, family = binomial(link = "logit"))$fitted.values
