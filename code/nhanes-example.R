# Demonstration of extrema fractional LP on NHANES data
# 
# devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))

library(tidyverse)
library(NHANES)
library(latex2exp)
library(decompsens)


## DAN AND QINGYUAN'S DATA ##
# load("../data/nhanes.fish.rda")
# The outcome of interest is log2(total blood mercury), measured in micrograms per liter; 
# the covariates include gender, age, income, whether income is missing and imputed, 
# race/ethnicity, education, smoking history, and the number of cigarettes smoked in the previous month.
# df <- nhanes.fish

## DATA FROM A RANDOM R PACKAGE ##
# data("mercury")
# # Guide: https://search.r-project.org/CRAN/refmans/submax/html/mercury.html
# df <- mercury


df <- NHANES

## Group and Outcome ##

group <- "Race1"
outcome <- "Depressed"
df %>% pull(!!sym(outcome)) %>% head()

df_cleaned <- df %>% 
  mutate(
    SexOrientation = case_when(
      SexOrientation %in% c("Bisexual", "Homosexual") ~ 1,
      SexOrientation == "Heterosexual" ~ 0
    ),
    Depressed = case_when(
      Depressed %in% c("Several", "Most") ~ 1,
      Depressed == "None" ~ 0
    ),
    Marijuana = case_when(
      Marijuana  == "Yes" ~ 1,
      Marijuana  == "No" ~ 0
    )
  ) %>% 
  mutate_at(vars(Race1), ~case_when(
    Race1 == "Mexican" | Race1 == "Hispanic" ~ "Latinx",
    .default = Race1
  ) %>% factor()) #%>% mutate_at(vars(!!sym(outcome)), ~as.numeric(.x) - 1)


df_cleaned %>% pull(!!sym(outcome)) %>% head()

df_cleaned %>% 
  group_by(!!sym(group)) %>% 
  drop_na(!!sym(group), !!sym(outcome)) %>%
  summarise(mean_out = mean(!!sym(outcome)),
            n = n())

print(paste0("We observe a disparity in ", outcome, " between ", group))

df_cleaned %>% pull(!!sym(outcome)) %>% table()

treatment <- "PhysActive"

df_cleaned %>% pull(!!sym(treatment)) %>% head()

df_cleaned <- df_cleaned %>% 
  mutate_at(vars(!!sym(treatment)), ~case_when(
    .x == "Yes" ~ 1,
    .x == "No" ~ 0
  ))

df_cleaned %>% pull(!!sym(treatment)) %>% head()

df_cleaned %>% 
  group_by(!!sym(group)) %>% 
  drop_na(!!sym(group), !!sym(outcome), !!sym(treatment)) %>%
  summarise(mean_trt = mean(!!sym(treatment)),
            n = n())

print(paste0("We observe a disparity in ", treatment, " between ", group))

df_cleaned %>% pull(!!sym(treatment)) %>% table()


## Covariates ##

names(df_cleaned)

allowable_covs <- c("Age", "Gender", "BMI")
non_allowable_covs <- c("Education", "HHIncome", "MaritalStatus", "Poverty")

df_analysis <- df_cleaned %>% 
  filter(!!sym(group) %in% c("Other", "Latinx")) %>% 
  mutate(Race1 = case_when(
    !!sym(group) == "Latinx" ~ 1,
    !!sym(group) == "Other" ~ 0
  )) %>% 
  dplyr::select(all_of(c(group, outcome, treatment, allowable_covs, non_allowable_covs))) %>% 
  drop_na()


XA <- model.matrix(~ ., data = df_analysis %>% dplyr::select(all_of(allowable_covs)))
XN <- model.matrix(~ ., data = df_analysis %>% dplyr::select(all_of(non_allowable_covs)))
X <- cbind(XA, XN[, -1])

Z <- df_analysis %>% dplyr::pull(treatment)
Y <- df_analysis %>% dplyr::pull(outcome)
G <- df_analysis %>% dplyr::pull(group)

w_rmpw = decompsens::estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN, trim = 0.05, allowable = FALSE)

w_rmpw[G == 1] %>% summary()

#w_rmpw <- pmin(w_rmpw, 5)

w_rmpw[G == 1] %>% summary()

# Counterfactual mean #
#mu_10 <- mean(w_rmpw * Y * G)
mu_10 <- sum(w_rmpw * Y * G) / sum(w_rmpw * G)
mu1 <- mean(Y[G == 1], na.rm = TRUE) 
mu0 <- mean(Y[G == 0], na.rm = TRUE)

mu1
mu_10
mu0

##################################################################
##                      Observed Disparity                      ##
##################################################################

true_disp <- mu1 - mu0
true_disp


##################################################################
##                      Residual Disparity                      ##
##################################################################
disp_residual <- lm(Y ~ G, weights = ifelse(G == 1, w_rmpw, 1))$coef[2]
disp_residual
mu_10 - mu0




#################################################################
##                     Disparity Reduction                     ##
#################################################################

df_reduction <- bind_rows(
  data.frame(Y, G, w_rmpw) %>% filter(G == 1) %>% mutate(D = 1),
  data.frame(Y, G, w_rmpw) %>% filter(G == 1) %>% mutate(D = 0)
)

disp_reduction <- lm(Y ~ D, data = df_reduction, weights = ifelse(D == 0, w_rmpw, 1))$coef[2]
disp_reduction
mu1 - mu_10
mu1/mu_10

mu_10/mu0

# Lam <- 1.12 # critical for residual disparity

Lam <- 1.12
estimand <- "res"


extrema <- getExtrema(G = G, Y = Y, 
                      gamma = log(Lam), w = w_rmpw, estimand = estimand, stab = TRUE)
extrema

bootci <- decompsens::boostrapCI(G, Z, Y, XA, XN, gamma = log(Lam), 
                                 alpha = 0.05, estimand = estimand,
                                 parallel = TRUE, B = 1000, stab = TRUE, 
                                 allow = FALSE, trim = 0.05)
bootci


### Amplification using Dan's Calibration ###

bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, Lambda = Lam)
str(bounds)

bounds[[1]] # muh - mu
bounds[[2]] # point estimate

maxbias <- max(abs(bounds[[1]]))
maxbias

####################################################
# Standardize observed covariates for control units
####################################################

X_G1 <- X[G == 1, -1]

# standardize X for control units
# center and make var = sd = 1
X_G1_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {scale(x)})

#################################################################
# Compute max coefficient among standardized observed covariates
#################################################################

coeffs_old <- numeric(ncol(X_G1))
for (var in 1:ncol(X_G1)) {
  coeffs_old[var] <- lm(Y[G==1] ~ X_G1[,-var] + X_G1_stnd[, var])$coef["X_G1_stnd[, var]"]
}
max(abs(coeffs_old))
# coeffs_old is equivalent to:

mod_matrix_y <- data.frame(y = Y[G==1], model.matrix(~ . - 1, data = X_G1_stnd %>% data.frame()))

coeffs <- lm(y ~ ., data = mod_matrix_y)$coef[-1]
max_betau_01 <- max(abs(coeffs))
max_betau_01


########################################################
# Compute maximum imbalance for standardized covariates
########################################################

X_stnd <- apply(X[, 2:ncol(X)], MARGIN = 2, FUN = function(x) {scale(x)})


imbal_stnd <- colMeans(X_stnd[G == 1, ]) - colMeans(X_stnd[G == 0, ])
max_imbal_stnd <- max(abs(imbal_stnd))
max_imbal_stnd
mean(abs(imbal_stnd))

####################################################
# Compute imbalance in covariates after weighting
####################################################


# imbal_stnd_weight <- colMeans(X_stnd[Z == 1,]) - colSums(X_stnd[Z == 0,]*w/sum(w))

####################################################
# Compute imbalance in covariates between G=0 and G=1_int after weighting for mu_10 
####################################################

w <- bounds[[3]]

imbal_stnd_weight <- colMeans(X_stnd[G == 0, ]) - colSums(X_stnd[G == 1, ] * (w[G == 1] / sum(w[G == 1]))) 
max_imbal_stnd_wt <- max(abs(imbal_stnd_weight))
max_imbal_stnd_wt
mean(abs(imbal_stnd_weight))

####################################################
# Compute imbalance in Group G = 1 after weighting for mu_10 
####################################################

w <- bounds[[3]]

imbal_stnd_weight <- colMeans(X_G1_stnd) - colSums(X_G1_stnd * (w[G == 1] / sum(w[G == 1])))
max_imbal_stnd_wt <- max(abs(imbal_stnd_weight))
max_imbal_stnd_wt
mean(abs(imbal_stnd_weight))



####################################################
# Get coordinates for strongest observed covariates to plot
####################################################
# get coefficients
coeff_df <- data.frame(
  covar = gsub("X_stnd[G == 1, ]", "", names(coeffs)),
  coeff = abs(as.numeric(coeffs)))

# get imbal
imbal_df <- data.frame(
  covar = names(imbal_stnd),
  imbal = abs(as.numeric(imbal_stnd)),
  imbal_wt = abs(as.numeric(imbal_stnd_weight)))
# merge coefficients and imbalance, arrange from largest to smallest by imbal * beta_u
strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>%
  dplyr::arrange(desc(coeff*imbal))

strongest_cov_df

#######
# Plot
#######

# function to compute beta_u from imbalance for plot
betauFun <- function(x) {
  maxbias/x
}


num_cov <- 4

x_orig = strongest_cov_df[1:num_cov,]$imbal_wt
y_orig = strongest_cov_df[1:num_cov,]$coeff

x <- c(x_orig,0, 0, max(x_orig))
y <- c(y_orig,0, max(y_orig),0)
df_plot <- matrix(c(x,y), ncol = 2)

hpts <- chull(df_plot)
hpts <- c(hpts, hpts[1])

X_cvx <- df_plot[hpts,]
mycurve1 <- as.data.frame(curve(betauFun, from=0.05, to=4)) # computes beta_u from imbalance

num_label <- 2

andy_ggplot_theme <- function(bs = 12) {
  theme_minimal(base_size = bs) + 
    # theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
    #       legend.title = element_blank(), legend.position = "none") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}



# Attempting to do this Melody-style

beta <- seq(0.05, 5.5, by = 0.01)
imbalance <- seq(0.001, 1, by = 0.005)
data.fit <- expand.grid(beta, imbalance)
names(data.fit) <- c("beta", "imbalance")

df_plot <- data.fit %>% 
  mutate(bias = beta * imbalance)


bins <- c(0.4, 0.8, 1.6)
text_bins <- bins[bins %% 1 == 0]

num_cov <- 5

p1 <- df_plot %>%
  ggplot(aes(x = imbalance, y = beta, z = bias)) +
  theme_bw() + 
  geom_contour(col="gray55", breaks = bins) + 
  metR::geom_text_contour(aes(z = bias), 
                          breaks = bins, 
                          stroke = 0.2, skip = 0) + 
  metR::geom_contour_fill(breaks = c(maxbias, 1000 * maxbias), fill='blue', alpha = 0.2) +
  geom_contour(breaks = c(maxbias), col='blue', size=1) 

p1

p1_full <- p1 + geom_point(data = strongest_cov_df[1:num_cov,], 
             aes(x = imbal, y = coeff, z = 0)) + 
  ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov,], 
                            aes(x = imbal, y = coeff, z = 0, label = covar),
                            nudge_y = 0.05, nudge_x = 0.5, label.padding = 0.1,
                            point.padding = 0.1)
p1_full

p1_full + 
  scale_x_continuous(name = "absolute standardized imbalance in U", limits = c(0, max(max_imbal_stnd, 1))) +
  scale_y_continuous(name = TeX("absolute $\\beta_u$"), limits = c(0, max(max_betau_01, betauFun(0.03))))

ggplot() +
  #geom_line(data = mycurve1, aes(x = x, y = y), color = "black", linewidth = 0.7) + # plots calibration curve
  stat_function(fun = betauFun, aes(), linewidth = 0.7) +
  geom_point(data = strongest_cov_df[1:num_cov,], 
             aes(x = imbal, y = coeff, colour = "unweighted observed covs")) + # plot unweighted imbalance points
  geom_point(data = strongest_cov_df[1:num_cov,], 
             aes(x = imbal_wt, y = coeff, colour = "weighted observed covs")) + # plot weighted imbalance points
  geom_hline(aes(yintercept=max(strongest_cov_df[1:num_cov,"coeff"]), colour = "unweighted observed covs"),
             linetype="dashed", alpha = 0.8) + # horizontal line at max beta_u
  geom_vline(aes(xintercept=max(strongest_cov_df[1:num_cov,"imbal"]), colour = "unweighted observed covs"),
             linetype="dashed", alpha = 0.8) + # vertical line at max imbalance
  geom_path(aes(x = X_cvx[,1], y = X_cvx[,2], colour ="weighted observed covs"), 
            alpha = 0.5) + 
  geom_polygon(aes(x = X_cvx[,1], y = X_cvx[,2]), 
               fill = "red", alpha = 0.2) +
  geom_text(data=strongest_cov_df[1:num_label,], 
            aes(x = imbal, y = coeff, colour = "unweighted observed covs",label=covar),
            size = 3, hjust = 0.1, vjust = 2.2) +
  #annotate("rect", xmin = 0, xmax = max(strongest_cov_df[1:num_cov,]$imbal_wt), 
  #          ymin = 0, ymax = max(strongest_cov_df[1:num_cov,]$coeff), alpha = .3) +
  scale_x_continuous(name = "absolute standardized imbalance in U", limits = c(0, max(max_imbal_stnd,3))) +
  scale_y_continuous(name = TeX("absolute $\\beta_u$"), limits = c(0, max(max_betau_01, betauFun(0.21)))) +
  #ggtitle(TeX(paste0("$\\beta_u$ vs. imbalance for  $\\Lambda$ = ", Lambda,", ", in_data, " data"))) +
  scale_color_manual(values = c("black","#00BFC4" , "#F8766D"), 
                      breaks = c("error", "unweighted observed covs", "weighted observed covs")) + 
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 10))
