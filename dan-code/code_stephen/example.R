# load package
source("bal_sens_att_ebal.R")
#install.packages("ebal")

library(parallel)

# load data
load("data/nhanes.fish.rda")

# treatment vector
Z <- as.numeric(nhanes.fish$fish.level == "high")
# design matrix
X_fish <- nhanes.fish[, c("gender", "age", "income", "income.missing",
                          "race", "education", "smoking.ever", "smoking.now")]
X_fish$race <- factor(X_fish$race)
X <- model.matrix(~ . - 1, X_fish)
# remove collinearity
X <- X[,colnames(X) != "race7"]
# outcome vector
Y <- log2(nhanes.fish$o.LBXTHG)

# Recreate Table 2 from Zhao et al. (2019) ATT SIPW with ebal

# Interval of point estimates and 90% confidence intervals for
# varying values of sensitivity parameter Lambda

# Lambda = 1 (no unmeasured confounding)
extrema.os(Z,X,Y, Lambda = 1, estimand = 'att')
bootsens.os(Z,X,Y, Lambda = 1, estimand = 'att', alpha = 0.1)

# Lambda = 1.65
extrema.os(Z,X,Y, Lambda = 1.65, estimand = 'att')
bootsens.os(Z,X,Y, Lambda = 1.65, estimand = 'att', alpha = 0.1)

# Lambda = 2.72
extrema.os(Z,X,Y, Lambda = 2.72, estimand = 'att')
bootsens.os(Z,X,Y, Lambda = 2.72, estimand = 'att', alpha = 0.1)

# Lambda = 7.39
extrema.os(Z,X,Y, Lambda = 7.39, estimand = 'att')
bootsens.os(Z,X,Y, Lambda = 7.39, estimand = 'att', alpha = 0.1)

# Lambda = 20.09
extrema.os(Z,X,Y, Lambda = 20.09, estimand = 'att')
bootsens.os(Z,X,Y, Lambda = 20.09, estimand = 'att', alpha = 0.1)
