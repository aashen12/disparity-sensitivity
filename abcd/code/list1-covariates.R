# Analysis of list 1

library(tidyverse)
library(decompsens)
library(parallel)
library(doParallel)
library(boot)
library(latex2exp)


set.seed(122357)

numCores <- parallel::detectCores()

doParallel::registerDoParallel(numCores)

df_x <- read_csv("../data/list1_X.csv")

# Allowable covariates 
options(na.action='na.pass')


#mediators <- c("peer_victimization")
mediators <- c("")
allowable_covs <- c("age", "sex", "sib_num", "sib_order")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

XA <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs)))
XN <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(non_allowable_covs)))


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

allow <- TRUE
trim <- 0

w <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=XA, XN=XN, trim = trim, allowable = allow)
summary(w)


X <- cbind(XA, XN)

# assess balance in X

lovePlot <- function(pre_weight, post_weight, title = "Covariate Balance") {
  data.frame(
    covariate = names(pre_weight),
    pre_weight = pre_weight,
    post_weight = post_weight
  ) %>% 
    pivot_longer(cols = c(pre_weight, post_weight), names_to = "balance", values_to = "mean_difference") %>% 
    ggplot(aes(x = mean_difference, y = covariate, color = balance)) +
    geom_point(size = 3) +
    theme_minimal(base_size = 16) +
    scale_color_manual(values = c("pre_weight" = "red", "post_weight" = "blue")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(title = title, x = "Mean Difference", y = "", color = "Sample") +
    scale_x_continuous(limits = c(-1, 1)) + # Set x-axis limits
    theme(legend.position = "bottom") + 
    theme(plot.title = element_text(hjust = 0.5))
}

X_stnd <- apply(X, 2, scale)
Xw <- apply(X_stnd, MARGIN = 2, FUN = function(x) {x * w / sum(w)})

pre_weight_Z <- colMeans(X_stnd[Z == 1, ]) - colMeans(X_stnd[Z == 0, ])
max(abs(pre_weight_Z))

post_weight_Z <- colMeans(Xw[Z == 1, ]) - colMeans(Xw[Z == 0, ])
max(abs(post_weight_Z))



loveZ <- lovePlot(pre_weight_Z, post_weight_Z, title = "Covariate Balance wrt Parental Support")
loveZ


pre_weight_G <- colMeans(X_stnd[G == 1, ]) - colMeans(X_stnd[G == 0, ])
max(abs(pre_weight_G))

post_weight_G <- colMeans(Xw[G == 1, ]) - colMeans(Xw[G == 0, ])
max(abs(post_weight_G))



loveG <- lovePlot(pre_weight_G, post_weight_G, title = "Covariate Balance wrt Sexual Minority Status")
loveG


# Assess balance in X for the sexual minority group

pre_weight_G1 <- colMeans(X_stnd[G == 1 & Z == 1, ]) - colMeans(X_stnd[G == 1 & Z == 0, ])
max(abs(pre_weight_G1))

post_weight_G1 <- colMeans(Xw[G == 1 & Z == 1, ]) - colMeans(Xw[G == 1 & Z == 0, ])
max(abs(post_weight_G1))

loveG1 <- lovePlot(pre_weight_G1, post_weight_G1, title = "Covariate Balance wrt Parental Support for Sex. Minorities")
loveG1

X_G1 <- X[G == 1, ] # [, -1] for intercept

# standardize X for group G = 1
X_G1_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {scale(x)})
ZG1 <- Z[G == 1]

pre_weight_GZ1 <- colMeans(X_G1_stnd[ZG1 == 1, ]) - colMeans(X_G1_stnd[ZG1 == 0, ])

wg1 <- w[G == 1]
X_G1_w_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {x * wg1 / sum(wg1)})

post_weight_GZ1 <- colMeans(X_G1_w_stnd[ZG1 == 1, ]) - colMeans(X_G1_w_stnd[ZG1 == 0, ])

loveGZ1 <- lovePlot(pre_weight_GZ1, post_weight_GZ1, title = "Covariate Balance wrt Parental Support for sex minorities")
loveGZ1







