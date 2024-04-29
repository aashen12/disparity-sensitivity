# Analysis of list 1
rm(list = ls())
# devtools::install_github("https://github.com/aashen12/decompsens")
# devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))

# this script computes the weights when one covariate is omitted at a time

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

outcome <- "attempt" # ideation or attempt

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

X_df <- cbind(df_allowable, df_non_allowable)

# Compute leave-one-out weights

covariate_names <- c(allowable_covs, non_allowable_covs, "sexF", "sexM")
covariate_names <- covariate_names[covariate_names != "sex"]

loco_weights <- map(
  seq_along(covariate_names),
  function(i) {
    cov <- covariate_names[i]
    if (cov %in% c(allowable_covs, "sexF", "sexM")) {
      contains_cov <- str_detect(names(df_allowable), cov)
      print(paste0("Covariate ", cov, " is allowable."))
      df_a_temp <- df_allowable[, !contains_cov]
      df_n_temp <- df_non_allowable
    } else {
      contains_cov <- str_detect(names(df_non_allowable), cov)
      print(paste0("Covariate ", cov, " is non-allowable."))
      df_n_temp <- df_non_allowable[, !contains_cov]
      df_a_temp <- df_allowable
    }
    weight_object <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=df_a_temp, XN=df_n_temp,
                                              trim = switch(outcome, "ideation" = 0.01, "attempt" = 0.05), 
                                              allowable = TRUE)
    w <- weight_object$w_rmpw
    e1 <- weight_object$e1
    e0 <- weight_object$e0
    list(w = w, e1 = e1, e0 = e0)
  }
); names(loco_weights) <- covariate_names


saveRDS(loco_weights, paste0("../data/loco_weights_", Z_method, "_", outcome, ".rds"))
print(paste0("Saved loco weights to ", paste0("../data/loco_weights_", Z_method, "_", outcome, ".rds")))



