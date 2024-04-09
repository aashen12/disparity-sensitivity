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

df_x <- read_csv("../data/list1_X.csv")

# Allowable covariates 
options(na.action='na.pass')


#mediators <- c("peer_victimization")
mediators <- c("")
allowable_covs <- c("age", "sex")
allowable_covs <- c("age", "sex", "sib_num", "sib_order")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)



df <- df_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute()
XA <- create_model_matrix(df_allowable)

df <- df_non_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute()
XN <- create_model_matrix(df_non_allowable)


# XA <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
#   data.frame() %>%
#   NAImpute()
# XA <- XA[, !grepl(":.*NA$", colnames(XA))
# XN <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% NAImpute()
# XN <- XN[, !grepl(":.*NA$", colnames(XN))]

X <- cbind(XA, XN)



dim(XA); dim(XN)

# Non-allowable covariates


df_yz <- read_csv("../data/list1_YZG.csv")

G <- df_yz$sex_min
Z <- df_yz$parent_accept

df_yz %>% group_by(sex_min) %>% 
  summarise(
    n = n(),
    mean_parent_accept = mean(parent_accept, na.rm = FALSE),
    sd_parent_accept = sd(parent_accept, na.rm = FALSE)
  )

allowable <- TRUE
trim <- 0.01

X <- cbind(XA, XN)
if (allowable) {
  e0 <- glm(Z ~ XA, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
  e1 <- glm(Z ~ X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
} else {
  e0 <- glm(Z ~ X, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
  e1 <- glm(Z ~ X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
}

e0 <- pmax(pmin(e0, 1 - trim), trim)
e1 <- pmax(pmin(e1, 1 - trim), trim)

#wr <- glm(G ~ 1, family = binomial, na.action = na.exclude)$fitted.values / glm(G ~ XA[, "age"] + XA[, "sexF"] + XA[, "sexM"], family = binomial, na.action = na.exclude)$fitted.values
wr <- 1

w1 = e0 / e1 * wr
w0 = (1 - e0) / (1 - e1) * wr
w = w1 * Z + w0 * (1 - Z)


summary(e0)
summary(e1)

hist(w1)
hist(w0)


# assess balance in X

lovePlot <- function(pre_weight, post_weight, num_covs = 100, title = "Covariate Balance") {
  subset_covs <- 2 * min(num_covs, length(pre_weight))
  pre_weight <- abs(pre_weight)
  post_weight <- abs(post_weight)
  if (any(names(pre_weight) != names(post_weight))) {
    stop("Pre and post weights must have the same names")
  }
  data.frame(
    covariate = names(pre_weight),
    pre_weight = pre_weight,
    post_weight = post_weight
  ) %>% 
    arrange(desc(pre_weight)) %>%
    mutate(covariate = fct_reorder(covariate, pre_weight)) %>% 
    pivot_longer(cols = c(pre_weight, post_weight), names_to = "balance", values_to = "mean_difference") %>% 
    mutate(allowability = ifelse(covariate %in% c(allowable_covs, "sexM", "sexF"), 1, 0)) %>% 
    mutate(covariate = fct_reorder(covariate, allowability)) %>% 
    slice(1:subset_covs) %>% 
    ggplot(aes(x = mean_difference, y = covariate, color = balance)) +
    geom_point(size = 6) +
    theme_minimal(base_size = 22) +
    scale_color_manual(values = c("pre_weight" = "red", "post_weight" = "blue")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(title = title, x = "Absolute Mean Difference", y = "", color = "Sample") +
    scale_x_continuous(limits = c(0, 1)) + # Set x-axis limits
    theme(legend.position = "bottom") + 
    theme(plot.title = element_text(hjust = 0.5)) #+ coord_flip()
}

X_plot <- X_plot <- cbind(model.matrix(~ . -1, data = df_x %>% select(all_of(allowable_covs))),
                          model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs)))) %>% NAImpute()
X_stnd <- apply(X_plot, 2, scale)



pre_weight_G <- colMeans(X_stnd[G == 1, ]) - colMeans(X_stnd[G == 0, ])
max(abs(pre_weight_G))

post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - 
  colSums(X_stnd[G == 0, ] * w[G == 0])  / sum(w[G == 0])
lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")

post_weight_G <- colSums(X_stnd[G == 1, ] * w[G == 1]) / sum(w[G == 1]) - colMeans(X_stnd[G == 0, ])
lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")


pre_weight_Z <- colMeans(X_stnd[Z == 1, ]) - colMeans(X_stnd[Z == 0, ])
max(abs(pre_weight_Z))

post_weight_Z <- colSums(X_stnd[Z == 1, ] * w[Z == 1] / sum(w[Z == 1])) - colSums(X_stnd[Z == 0, ] * w[Z == 0] / sum(w[Z == 0]))
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





