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





