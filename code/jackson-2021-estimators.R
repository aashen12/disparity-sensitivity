set.seed(12)

library(tidyverse)
library(bannerCommenter)


# Number of observations
n <- 10

# Generating data
R <- rbinom(n, 1, 0.5)
Ay1 <- rnorm(n)  # Background covariate 1
Ay2 <- rnorm(n)  # Background covariate 2

# generate M with values dependent on R

# category[race == 1] <- sample(c("A", "B", "C"), sum(race == 1), replace = TRUE, prob = c(0.5, 0.3, 0.2))
# category[race == 0] <- sample(c("A", "B", "C"), sum(race == 0), replace = TRUE, prob = c(0.2, 0.3, 0.5))

M <- character(n)
M[R == 1] <- sample(c("L1", "L2", "L3"), sum(R == 1), replace = TRUE, prob = c(0.1, 0.2, 0.7))
M[R == 0] <- sample(c("L1", "L2", "L3"), sum(R == 0), replace = TRUE, prob = c(0.2, 0.3, 0.5))

cat_coefs <- ifelse(M == "L1", 2, 
                    ifelse(M == "L2", 4, 7))

Y <- 2 + 7 * R + 0.2 * Ay1 + 0.5 * Ay2 + cat_coefs + rnorm(n)

df <- tibble(Y, R, Ay1, Ay2, M)

pscore <- glm(R ~ Ay1 + Ay2, family = "binomial", data = df)$fitted.values

df_all <- df %>% select(Y, R, Ay1, Ay2, M) %>% 
  arrange(desc(R)) %>% 
  mutate(D = ifelse(R == 1, 1, 0))
# if race == 1, D == 1 for now since we have not duplicated R == 1 yet

# create w_1 and w_0
df_all <- df_all %>% mutate(w_1 = mean(R == 1) / pscore, w_0 = mean(R == 0) / (1 - pscore))

# create w_r
df_all <- df_all %>% mutate(w_r = ifelse(R == 1, w_1, w_0))

# create w_rmpw
df_all <- df_all %>% mutate(w_rmpw = w_1 * mean(M[R == 0] == "L3") / mean(M[R == 1] == "L3"))

# duplicate rows with R == 1 and append to existing df_all
df_all <- df_all %>% 
  filter(R == 1) %>% 
  mutate(D = 0) %>%
  # add stacked R == 1 group with D == 0
  bind_rows(df_all) %>% 
  mutate(w_reduction = ifelse(D == 1, w_r, w_rmpw)) %>% 
  # if D == 1, they are in the non-intervened black group so they are weighted by the usual weight
  # if D == 0, they are in the intervened black group so they are weighted by the rmpw weight
  mutate(w_residual = ifelse(R == 0, w_r, w_rmpw))
# if R == 0, they are white so they are weighted by the usual weight
# if R == 1, they are intervened on so they are weighted by the rmpw weight


# compute each estimator
mu1hat <- with(df_all, sum(w_r * R * D * Y) / sum(w_r * D * R))
mu0hat <- with(df_all, sum(w_r * (1 - R) * Y) / sum(w_r * (1 - R)))
mu_10 <- with(df_all, sum(w_rmpw * R * (1 - D) * Y) / sum(w_rmpw * R * (1 - D)))


# compute each estimator: rewrite based on spi feedback
mu1hat2 <- with(df_all %>% filter(D == 0), sum(w_r * R * Y) / sum(w_r * R))
mu0hat2 <- with(df_all %>% filter(D == 0), sum(w_r * (1 - R) * Y) / sum(w_r * (1 - R)))
mu_102 <- with(df_all %>% filter(D == 0), sum(w_rmpw * R * Y) / sum(w_rmpw * R))

# observed disparity
mu1hat - mu0hat
mu1hat2 - mu0hat2

mod_obs <- lm(Y ~ R, data = df_all %>% filter(D == 0), weights = w_r)
mod_obs$coefficients["R"]



# disparity reduction
mu1hat - mu_10
mu1hat2 - mu_102

mod_reduction <- lm(Y ~ D, data = df_all %>% filter(R == 1), weights = w_reduction)
mod_reduction$coefficients["D"]



# residual disparity
mu_10 - mu0hat
mu_102 - mu0hat2

mod_resid <- lm(Y ~ R, data = df_all %>% filter(D == 0), weights = w_residual)
mod_resid$coefficients["R"]



