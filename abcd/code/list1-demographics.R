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
library(table1)
library(lsr)
library(knitr)
library(printr)

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

outcome <- "ideation" # ideation or attempt


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
# allowable_covs <- c("age", "sex", "sib_order", "sib_num")
allowable_covs <- c("age", "sex")
#allowable_covs <- c("age", "sex", "family_conflict")
#allowable_covs <- c("age", "sex", "sib_num", "sib_order", "income", "adi")
non_allowable_covs <- setdiff(names(df_x)[!names(df_x) %in% mediators], allowable_covs)

df_allowable <- model.matrix(~ .^2 -1, data = df_x %>% select(all_of(allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

df_non_allowable <- model.matrix(~ . -1, data = df_x %>% select(all_of(non_allowable_covs))) %>% 
  data.frame() %>% NAImpute() %>% tibble()

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (%s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f%%)", FREQ, PCT))))
}

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

demog_df <- full_join(df_x, df_yz, by = "src_subject_id")
demog_df$sex_min <- factor(demog_df$sex_min, levels = c(0, 1), labels = c("Heterosexual", "Sexual Minority"))
demog_df$parent_accept <- factor(demog_df$parent_accept, levels = c(0, 1), labels = c("Poor Parental Acceptance", "Superior Parental Acceptance"))
demog_df$sex <- factor(demog_df$sex, levels = c("F", "M"), labels = c("Female", "Male"))
demog_df$family_mental_health <- as.logical(demog_df$family_mental_health)
demog_df$income <- factor(demog_df$income, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
label(demog_df$parent_accept) <- "Parental Acceptance"
label(demog_df$sex_min) <- "Sexual Minority Status"
label(demog_df$age) <- "Age"
units(demog_df$age) <- "years"
label(demog_df$sex) <- "Sex Assigned at Birth"
label(demog_df$sib_num) <- "Number of Siblings"
label(demog_df$sib_order) <- "Birth Order"
label(demog_df$income) <- "Household Income"
units(demog_df$income) <- "scale of 1-10, 1 is lowest income bracket"
label(demog_df$adi) <- "Area Deprivation Index"
units(demog_df$adi) <- "numeric measurement, larger is more deprived"
label(demog_df$family_conflict) <- "Family Conflict"
units(demog_df$family_conflict) <- "numeric measurement, larger is more conflict"
label(demog_df$family_mental_health) <- "Previous Suicide Attempt in Family"
label(demog_df$peer_victimization.x) <- "Peer Victimization"
units(demog_df$peer_victimization.x) <- "numeric measurement, larger is more severe"
label(demog_df$school_safety) <- "School Safety"
units(demog_df$school_safety) <- "numeric measurement, larger is safer"
label(demog_df$neighborhood_safety) <- "Neighborhood Safety"
units(demog_df$neighborhood_safety) <- "numeric measurement, larger is safer"
label(demog_df$structural_stigma) <- "Structural Stigma"
units(demog_df$structural_stigma) <- "State level indicators of sexism from survey and implicit bias measures"

t1 <- table1(~age + sex + sib_num + sib_order + income + adi + 
               family_conflict + family_mental_health + peer_victimization.x + 
               school_safety + neighborhood_safety + adi + structural_stigma | sex_min*parent_accept, 
             data = demog_df, 
             render.continuous=my.render.cont, 
             render.categorical=my.render.cat)
x <- t1

x
