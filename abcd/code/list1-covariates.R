# Analysis of list 1
rm(list = ls())

library(tidyverse)
library(decompsens)
library(parallel)
library(doParallel)
library(boot)
library(cobalt)
library(latex2exp)

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
mediators <- c("src_subject_id")
allowable_covs <- c("age", "sex", "sib_num", "sib_order")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

df_non_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
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

allow <- TRUE

weight_object <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=df_allowable, XN=df_non_allowable,
                              trim = switch(outcome, "ideation" = 0.01, "attempt" = 0.05), 
                              allowable = allow)


w <- weight_object$w_rmpw
e1 <- weight_object$e1
e0 <- weight_object$e0

w1 <- glm(G ~ 1, data = df_allowable, family = binomial(link = "logit"))$fitted.values / 
  glm(G ~ ., data = df_allowable, family = binomial(link = "logit"))$fitted.values
w0 <- glm(1 - G ~ 1, data = df_allowable, family = binomial(link = "logit"))$fitted.values / 
  glm(1 - G ~ ., data = df_allowable, family = binomial(link = "logit"))$fitted.values

summary(w)

df_yz <- df_yz %>% mutate(w_rmpw = w,
                          e1 = e1,
                          e0 = e0,)
write_csv(df_yz, paste0("../data/list1_YZGW_", Z_method, "_", outcome, ".csv"))

message(paste0("CSV file for Z method ", Z_method, " and outcome ", outcome, " has been written."))



X_plot <- cbind(model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))),
                model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs)))) %>% NAImpute()

# X_plot <- cbind(model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))),
#                 model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs)))) %>% NAImpute()


X_stnd <- apply(X_plot, 2, scale) %>% data.frame()


unweighted_diffs <- lapply(1:ncol(X_stnd), function(i) {
  tapply(X_stnd[, i], list(G, Z), mean)
})


weighted_diffs <- lapply(1:ncol(X_stnd), function(i) {
  tapply(X_stnd[, i] * w, list(G, Z), mean)
})


df_balance <- data.frame(
  X_stnd, G, Z, w
)

df_balance %>% 
  group_by(G, Z) %>% 
  summarise(across(colnames(X_stnd), mean))




pre_weight_G <- colMeans(X_stnd[G == 1, ]) - colMeans(X_stnd[G == 0, ])
max(abs(pre_weight_G))

# post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - 
#   colSums(X_stnd[G == 0, ] * w[G == 0])  / sum(w[G == 0])
# lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")

post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - colMeans(X_stnd[G == 1, ])
lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance between SM and weighted SM")


pre_weight_Z <- colMeans(X_stnd[Z == 1, ]) - colMeans(X_stnd[Z == 0, ])
pre_weight_Z <- colSums(X_stnd[Z == 1, ] * w1[Z == 1] / sum(w1[Z == 1])) - 
  colSums(X_stnd[Z == 0, ] * w0[Z == 0] / sum(w0[Z == 0]))
max(abs(pre_weight_Z))

post_weight_Z <- colSums(X_stnd[Z == 1, ] * w[Z == 1] / sum(w[Z == 1])) - colMeans(X_stnd[Z == 0, ])
max(abs(post_weight_Z))

#XG1 <- X_stnd[G == 1, ]
XG1_stnd <- apply(X_plot[G == 1, ], 2, scale)
XG1_w <- apply(XG1_stnd, 2, function(x) x * w[G == 1] / sum(w[G == 1]))

# compute scaling factors 
sf <- ((e0 - e1) / (1 - e1))[G == 1] #%>% abs()
length(sf)
dim(XG1_stnd)
X_sf <- apply(XG1_stnd, 2, function(x) x * sf / sum(abs(sf)))

ZG1 <- Z[G == 1]
sum(ZG1)

sf_Z1 <- ((e0 - e1) / (e1 - e1^2))[G == 1 & Z == 1] #%>% abs()
dim(XG1_stnd[ZG1 == 1, ])
length(sf_Z1)
X_sf_Z1 <- apply(XG1_stnd[ZG1 == 1, ], 2, function(x) x * sf_Z1 / sum(abs(sf_Z1)))


pre_weight_Z <- colMeans(XG1_stnd) - colMeans(XG1_stnd[ZG1 == 1, ])
post_weight_Z <- colSums(XG1_w) - colSums(XG1_w[ZG1 == 1, ])
#post_weight_Z <- colMeans(XG1_stnd) - colSums(XG1_w[ZG1 == 1, ])
post_weight_Z <- colSums(X_sf) - colSums(X_sf_Z1)

loveZ <- lovePlot(pre_weight_Z, post_weight_Z, 
                  title = "Covariate Balance between all SM and treated SM")
loveZ

