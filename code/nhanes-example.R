# Demonstration of extrema fractional LP on NHANES data

library(tidyverse)
library(submax)

source("extrema_and_bootstrap_demo.R")

## DAN AND QINGYUAN'S DATA ##
# load("../data/nhanes.fish.rda")
# The outcome of interest is log2(total blood mercury), measured in micrograms per liter; 
# the covariates include gender, age, income, whether income is missing and imputed, 
# race/ethnicity, education, smoking history, and the number of cigarettes smoked in the previous month.
# df <- nhanes.fish

## DATA FROM A RANDOM R PACKAGE ##
data("mercury")
# Guide: https://search.r-project.org/CRAN/refmans/submax/html/mercury.html
df <- mercury

head(df)

names(df)
unique(df$black)

G <- df$black
Y <- df$methylmercury
Z <- df$fish
X <- df$female

e0 <- glm(fish ~ female, data = df, family = binomial, weights = 1 - black)$fitted.values
e1 <- glm(fish ~ female, data = df, family = binomial, weights = black)$fitted.values

w_rmpw <- e1 * Z + e0 * (1 - Z)

mu1 <- mean(Y[G == 1])
mu0 <- mean(Y[G == 0])
mu_10 <- sum(w_rmpw * Y * G) / sum(w_rmpw * G)

mu1
mu0
mu_10

mu1 - mu_10
mu_10 - mu0





estimand <- "res"
Lam <- 1
alpha <- 0.05

extrema <- getExtrema(G, Y, gamma = log(Lam), w = w_rmpw, estimand = estimand)
extrema

confint <- boostrapCI(G, Y, gamma = log(Lam), w = w_rmpw, 
                      alpha = alpha, estimand = estimand, B = 8000)
confint

