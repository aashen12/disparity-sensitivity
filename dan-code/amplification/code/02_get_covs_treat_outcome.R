library(tidyverse)
library(TukeySens)
library(foreign)
##############
# load inputs
##############
input_df <- readRDS("../inputs/input_df.rds")

############
# Read data
############
if (input_df$data == 'rhc') {
  rhc <- read.csv("../rhc/data/rhc.csv")
} else if (input_df$data %in% c('lalonde', 'lalonde interactions')) {
  # read data
  # https://users.nber.org/~rdehejia/data/.nswdata2.html
  nsw <- read.dta("../lalonde/data/nsw_dw.dta")
  cps_control <- read.dta("../lalonde/data/cps_controls.dta")
} else if (input_df$data == 'school') {
  
} else if (input_df$data == 'fish') {
  # https://cran.r-project.org/web/packages/CrossScreening/index.html
  # https://cran.r-project.org/src/contrib/Archive/CrossScreening/
  load("../fish/data/nhanes.fish.rda")
}

########################################
# Define treatment, covariates, outcome
########################################

if (input_df$data == 'rhc') {
  # Treatment: "use of right heart catheterization (RHC) during first 24 hours of care in
  #            the intensive care unit (ICU)"
  # treatment vector
  Z <- as.numeric(rhc$swang1 == "RHC")
  
  # Outcome: "patient survival time, hospital and ICU length of stay, hospital costs, 
  #           and intensity of care"
  
  # outcome vector
  Y <- as.numeric(rhc$dth30 == "Yes")
  
  # Covariates: "age, sex, race (black, white, other), years of education, income, 
  #   type of medical insurance (private, Medicare, Medicaid, private and Medicare, 
  # Medicare and Medicaid, or none), primary disease category, secondary disease category,
  # 12 categories of admission diagnosis,
  # ADL and DASI 2 weeks before admission,
  # do-not-resuscitate status on day 1,
  # cancer (none, localized, metastatic),
  # SUPPORT model estimate of the probability of surviving 2 months, 
  # 17 acute physiology component of the AP ACHE III score, 16 Glasgow Coma Score, weight,
  # temperature, mean blood pressure, respiratory
  # rate, heart rate, PaO:/FI02 ratio,
  # PaC02, pH, WBC count, hematocrit,
  # sodium, potassium, creatinine, bilirubin,
  # albumin, urine output, and 13 categories
  # of comorbid illness"
  
  # and 13 categories of comorbid illness
  # ADL and DASI 2 weeks before admission,
  # do-not-resuscitate status on day 1,
  # 16 Glasgow Coma Score
  
  
  # design matrix
  # X_df <- rhc[, c("age", "sex","race", "edu", "income", "ninsclas", "cat1", "cat2",
  #                 "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx",
  #                 "gibledhx", "malighx", "immunhx", "transhx", "amihx", "adld3p","das2d3pc","dnr1", "ca", "surv2md1", "aps1", 
  #                 "scoma1","wtkilo1", "temp1", "meanbp1", "resp1", "hrt1", "pafi1", "paco21","ph1","wblc1", "hema1",
  #                 "sod1", "pot1", "crea1", "bili1", "alb1", "urin1",
  #                 "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho" )]
  
  # remove columns with NAs
  # X_df <- rhc[, c("age", "sex","race", "edu", "income", "ninsclas", "cat1",
  #                 "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx",
  #                 "gibledhx", "malighx", "immunhx", "transhx", "amihx","das2d3pc","dnr1", "ca", "surv2md1", "aps1", 
  #                 "scoma1","wtkilo1", "temp1", "meanbp1", "resp1", "hrt1", "pafi1", "paco21","ph1","wblc1", "hema1",
  #                 "sod1", "pot1", "crea1", "bili1", "alb1",
  #                 "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho" )]
  
  # remove columns with NAs
  # get rid of trauma and ortho. too many zeros so bootstrap samples sometimes have all 0s
  X_df <- rhc[, c("age", "sex","race", "edu", "income", "ninsclas", "cat1",
                  "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx", "liverhx",
                  "gibledhx", "malighx", "immunhx", "transhx", "amihx","das2d3pc","dnr1", "ca", "surv2md1", "aps1", 
                  "scoma1","wtkilo1", "temp1", "meanbp1", "resp1", "hrt1", "pafi1", "paco21","ph1","wblc1", "hema1",
                  "sod1", "pot1", "crea1", "bili1", "alb1",
                  "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps" )]
  
  
  X <- model.matrix(~ . - 1, X_df)
  X <- X[,!grepl("sexMale",colnames(X))]
  
} else if (input_df$data == 'nhanes') {
  # treatment vector
  Z <- NHANES$trt_dbp
  # design matrix
  X_df <- NHANES[,!names(NHANES) %in% c("trt_dbp", "ave_dbp")]
  X <- model.matrix(~ . - 1, X_df)
  # outcome vector
  Y <- NHANES$ave_dbp
} else if (input_df$data == 'lalonde') {
  # nsw treatment group
  nsw_treat <- filter(nsw, treat == 1)
  rm(nsw)
  
  # combine
  nsw_combined <- rbind(nsw_treat, cps_control)
  
  # create employment status var
  
  nsw_combined$unemployed_74 <- as.numeric(nsw_combined$re74 == 0)
  nsw_combined$unemployed_75 <- as.numeric(nsw_combined$re75 == 0)
  
  # create treatment and outcome vectors, design matrix
  
  # treatment vector
  Z <- nsw_combined$treat
  # design matrix
  X_nsw <- dplyr::select(nsw_combined, -c(treat, data_id, re78))
  X <- model.matrix(~ . - 1, X_nsw)
  # outcome vector
  Y <- nsw_combined$re78
} else if (input_df$data == 'lalonde interactions') {
  
  # nsw treatment group
  nsw_treat <- filter(nsw, treat == 1)
  #mean(nsw[nsw$treat==1,]$re78) - mean(nsw[nsw$treat==0,]$re78 )
  rm(nsw)
  
  # combine
  nsw_combined <- rbind(nsw_treat, cps_control)
  
  # create employment status var (earnings = 0 -> unemployed)
  nsw_combined$unemployed_74 <- as.numeric(nsw_combined$re74 == 0)
  nsw_combined$unemployed_75 <- as.numeric(nsw_combined$re75 == 0)
  
  # create treatment and outcome vectors, design matrix
  
  # treatment vector
  Z <- nsw_combined$treat
  # design matrix
  X_nsw <- dplyr::select(nsw_combined, -c(treat, data_id, re78))
  
  # From Hainmueller:
  # 10 raw variables, all their pairwise one-way interactions, as well as
  # squared terms for the continuous variables age and years of education. 
  # Overall, this results in 52 covariate combinations
  
  # all interactions
  X_nsw_all <- as.data.frame(model.matrix( ~.^2, data=X_nsw)[,-1])
  X_nsw_all <- dplyr::select(X_nsw_all, -c("black:hispanic", "re74:re75", "education:nodegree",
                                           "re74:unemployed_74", "re75:unemployed_75"))
  # squared terms for the continuous variables age and years of education
  X_nsw_all$age_sq <- X_nsw_all$age * X_nsw_all$age
  X_nsw_all$education_sq <- X_nsw_all$education * X_nsw_all$education
  X <- as.matrix(X_nsw_all)

  # outcome vector
  Y <- nsw_combined$re78
} else if (input_df$data == 'school') {
  
} else if (input_df$data == 'fish') {
  # treatment vector
  Z <- as.numeric(nhanes.fish$fish.level == "high")
  # design matrix
  X_fish <- nhanes.fish[, c("gender", "age", "income", "income.missing",
                       "race", "education", "smoking.ever", "smoking.now")]
  X_fish$race <- factor(X_fish$race)
  X <- model.matrix(~ . - 1, X_fish)
  # remove collinearity
  X <- X[,colnames(X) != "race7"]
  # convert gender into 0/1 isntead of 1/2
  X[,"gender"] <- X[,"gender"] - 1
  # outcome vector
  Y <- log2(nhanes.fish$o.LBXTHG)
}

# save
saveRDS(Z, file = paste0("../", input_df$data, "/intermediate/Z_", input_df$data, ".rds"))
saveRDS(Y, file = paste0("../", input_df$data, "/intermediate/Y_", input_df$data, ".rds"))
saveRDS(X, file = paste0("../", input_df$data, "/intermediate/X_", input_df$data, ".rds"))
