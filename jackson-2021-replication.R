# Setting the seed for reproducibility
set.seed(123)

library(tidyverse)
library(bannerCommenter)



# Number of observations
n <- 600

# Generating data
R <- sample(c(0, 1), n, replace = TRUE)
Ay1 <- rnorm(n)  # Background covariate 1
Ay2 <- rnorm(n)  # Background covariate 2
Y <- 2 + 7 * R + 0.2 * Ay1 + 0.5 * Ay2 + rnorm(n)

df <- tibble(Y, R, Ay1, Ay2)
df_orig <- df

banner("Estimator 1:", "Observed Disparity", emph = TRUE)

############################################################################
############################################################################
###                                                                      ###
###                             ESTIMATOR 1:                             ###
###                          OBSERVED DISPARITY                          ###
###                                                                      ###
############################################################################
############################################################################


pscore <- glm(R ~ Ay1 + Ay2, family = "binomial", data = df)$fitted.values 

glm(R ~ 1, family = "binomial", data = df)$fitted.values 

# Defining weights
w1 <- mean(df$R == 1) / pscore
w0 <- mean(df$R == 0) / (1 - pscore)
w_r <- ifelse(R == 1, w1, w0) 
df$w1 <- w1
df$w_r <- w_r


# observed disparity without weighting
mean(df$Y[df$R == 1]) - mean(df$Y[df$R == 0])



# observed disparity with weighting

# definitely correct

# mu1hat <- sum(df$w_r * df$R * df$Y) / sum(df$w_r * df$R)
# mu0hat <- sum(df$w_r * (1 - df$R) * df$Y) / sum(df$w_r * (1 - df$R))

mu1hat <- sum(df$w_r[df$R == 1] * df$Y[df$R == 1]) / sum(df$w_r[df$R == 1])
mu0hat <- sum(df$w_r[df$R == 0] * df$Y[df$R == 0]) / sum(df$w_r[df$R == 0])

tauhat_as <- mu1hat - mu0hat

# observed disparity from prediction
#mean(mod1) - mean(mod0)

# Fitting the model
model <- lm(Y ~ R, data = df, weights = w_r)
# model <- lm(Y ~ R - 1, weights = w_r) equiv to model2 <- lm(Yw ~ Rw - 1) 

# Displaying the results
model$coefficients

# "ground truth"
tauhat_jj <- model$coefficients["R"]

tauhat_as
tauhat_jj

#############################

banner("Estimator 2:", "Disparity Reduction", emph = TRUE)

###########################################################################
###########################################################################
###                                                                     ###
###                            ESTIMATOR 2:                             ###
###                         DISPARITY REDUCTION                         ###
###                                                                     ###
###########################################################################
###########################################################################

# \mu_1 - \mu_1^0

M <- sample(c("L1", "L2", "L3"), n, replace = TRUE) %>% factor()

df_reduc <- df %>% mutate(M = M, .after = "Ay2") %>% 
  mutate(w_rmpw = w1 * mean(M[R == 0] == "L3") / mean(M[R == 1] == "L3"))

df_1 <- df_G <- df_reduc %>% filter(R == 1)

df_1 <- df_1 %>% mutate(D = 1)
df_G <- df_G %>% mutate(D = 0)
df_1G <- bind_rows(df_1, df_G)

df_1G <- df_1G %>% mutate(w_reg = ifelse(D == 1, w_r, w_rmpw))

# Fitting the model: ground truth
model <- lm(Y ~ D, data = df_1G, weights = w_reg)
tauhat_jj <- model$coefficients["D"]
tauhat_jj

mu1hat <- sum(df_1G$w_reg[df_1G$D == 1] * df_1G$Y[df_1G$D == 1]) / sum(df_1G$w_reg[df_1G$D == 1])
mu_10 <- sum(df_1G$w_reg[df_1G$D == 0] * df_1G$Y[df_1G$D == 0]) / sum(df_1G$w_reg[df_1G$D == 0])

mu1hat - mu_10

banner("Estimator 3:", "Residual Disparity", emph = TRUE)

############################################################################
############################################################################
###                                                                      ###
###                             ESTIMATOR 3:                             ###
###                          RESIDUAL DISPARITY                          ###
###                                                                      ###
############################################################################
############################################################################

df_resid <- df %>% mutate(M = M, .after = "Ay2") %>% 
  mutate(w_rmpw = w1 * mean(M[R == 0] == "L3") / mean(M[R == 1] == "L3")) %>% 
  mutate(w_reg = ifelse(R == 0, w_r, w_rmpw))
model <- lm(Y ~ R, data = df_resid, weights = w_reg)
tauhat_jj <- model$coefficients["R"]
tauhat_jj

mu_10 - mu0hat

banner("Putting it all together", emph = TRUE)

###########################################################################
###########################################################################
###                                                                     ###
###                       PUTTING IT ALL TOGETHER                       ###
###                                                                     ###
###########################################################################
###########################################################################








