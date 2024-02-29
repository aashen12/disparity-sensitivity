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
non_allowable_covs <- c("Education", "MaritalStatus", "Poverty")

df_analysis <- df_cleaned %>% 
  filter(!!sym(group) %in% c("White", "Latinx")) %>% 
  mutate(Race1 = case_when(
    !!sym(group) == "Latinx" ~ 1,
    !!sym(group) == "White" ~ 0
  )) %>% 
  dplyr::select(all_of(c(group, outcome, treatment, allowable_covs, non_allowable_covs))) %>% 
  mutate_at("Education",
            ~case_when(
              .x == "8th Grade" ~ 1,
              .x == "9 - 11th Grade" ~ 2,
              .x == "High School" ~ 3,
              .x == "Some College" ~ 4,
              .x == "College Grad" ~ 5
            )) %>% drop_na()


XA <- model.matrix(~ . -1, data = df_analysis %>% dplyr::select(all_of(allowable_covs)))
XN <- model.matrix(~ . -1, data = df_analysis %>% dplyr::select(all_of(non_allowable_covs)))
X <- cbind(XA, XN)

Z <- df_analysis %>% dplyr::pull(treatment)
Y <- df_analysis %>% dplyr::pull(outcome)
G <- df_analysis %>% dplyr::pull(group)

w_rmpw = decompsens::estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN, trim = 0.1, allowable = FALSE)

w_rmpw[G == 1] %>% summary()

#w_rmpw <- pmin(w_rmpw, 5)


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

mu1/mu0
mu1/mu_10
mu_10/mu0

# Lam <- 1.12 # critical for residual disparity



lam_resid <- 1.16 # critical for residual disparity alpha = 0.1
estimand <- "res"

w_rmpw <- decompsens::estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN, trim = 0.1, allowable = FALSE)

extrema <- getExtrema(G = G, Y = Y, 
                      gamma = log(lam_resid), w = w_rmpw, estimand = "po", stab = TRUE)
extrema

bootci <- decompsens::bootstrapCI(G, Z, Y, XA, XN, gamma = log(lam_resid), 
                                 alpha = 0.1, estimand = estimand,
                                 parallel = TRUE, B = 2000, stab = TRUE, 
                                 allowable = FALSE, trim = 0.1)
bootci

Lam <- lam_resid

bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, Lambda = Lam, trim = 0.1)

trim = 0.1
Lambda=Lam
bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, Lambda, 
                                    trim = trim, allowable = F, stab = T)
maxbias <- max(abs(bounds[[1]]))
message("Max bias: ", maxbias)

X_G1 <- X[G == 1, ]
X_G1_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {
  scale(x)
})
mod_matrix_y <- data.frame(y = Y[G == 1], model.matrix(~. - 
                                                         1, data = data.frame(X_G1_stnd)))
coeffs <- lm(y ~ ., data = mod_matrix_y)$coef[-1]
max_betau <- max(abs(coeffs), na.rm = TRUE)
ZG1 <- Z[G == 1]
imbal_stnd <- colMeans(X_G1_stnd[ZG1 == 1, ]) - colMeans(X_G1_stnd[ZG1 == 
                                                                     0, ])
max_imbal_stnd <- max(abs(imbal_stnd), na.rm = TRUE)
w <- bounds[[3]]
wg1 <- w[G == 1]
X_G1_w_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {
  x * wg1/sum(wg1)
})
imbal_stnd_weight <- colMeans(X_G1_w_stnd[ZG1 == 1, ]) - 
  colMeans(X_G1_w_stnd[ZG1 == 0, ])
max_imbal_stnd_wt <- max(abs(imbal_stnd_weight), na.rm = TRUE)
coeff_df <- data.frame(covar = gsub("X_stnd[G == 1, ]", "", 
                                    names(coeffs)), coeff = abs(as.numeric(coeffs)))
imbal_df <- data.frame(covar = names(imbal_stnd), imbal = abs(as.numeric(imbal_stnd)), 
                       imbal_wt = abs(as.numeric(imbal_stnd_weight)))
strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, 
                                      by = "covar") %>% dplyr::arrange(desc(coeff * imbal))


amplification <- decompsens::informalAmplify(G, Z, XA, XN, Y, Lambda = Lam, trim = 0.1)

strongest_cov_df <- amplification[[1]] %>% drop_na()
max_imbal_stnd <- amplification$max_imbal_stnd
max_betau_01 <- amplification$max_beta
str(amplification)

maxbias <- max(abs(bounds[[1]]))

max_betau_01 <- max(strongest_cov_df$coeff)
max_imbal <- max(strongest_cov_df$imbal)


round_down <- function(a) {
  floor(a * 100) / 100
}

beta <- seq(max(0, min(strongest_cov_df$coeff) - 0.2), max(strongest_cov_df$coeff) + 0.22, by = 0.005)
imbalance <- seq(0.001, max(strongest_cov_df$imbal) + 0.1, by = 0.005)
data.fit <- expand.grid(beta, imbalance)


names(data.fit) <- c("beta", "imbalance")

df_plot <- data.fit %>% 
  mutate(bias = abs(beta * imbalance))

summary(df_plot$bias)
maxbias

bins <- seq(max(df_plot$bias) - 1.6 * sd(df_plot$bias),  max(df_plot$bias) + 1 * sd(df_plot$bias), length.out = 4) %>% 
  round_down()

num_cov <- 8 # specify 2x what you want
psize <- 6.5

p1 <- df_plot %>%
  ggplot(aes(x = imbalance, y = beta, z = bias)) +
  theme_bw() + 
  geom_contour(col="gray55", breaks = bins) + 
  metR::geom_text_contour(aes(z = bias), 
                          stroke = 0.2, skip = 0, breaks = bins) + 
  metR::geom_contour_fill(breaks = c(maxbias, 1000 * maxbias), fill='powderblue', alpha = 0.5) +
  geom_contour(breaks = c(maxbias), col='blue', linewidth = 1) + 
  metR::geom_text_contour(aes(z = maxbias), stroke = 0.2) + 
  geom_point(x = 0.3, y = maxbias / 0.3, size = psize - 0.2, color = "black") +
  geom_text(x = 0.3 + 0.033 , y = maxbias / 0.3 + 0.01, 
            label = paste0("Bias: ", round(maxbias, 3)), size = psize, color = "black")

# geom_contour(data = data.frame(x = seq(0.01, 0.25, by = 0.01), y = seq(0.01, 2, by = 0.01)) %>% 
#                mutate(z = x * y),
#              aes(x = x, y = y, z = z), breaks = c(maxbias), col='blue', linewidth = 1, inherit.aes = F) +

p1

strongest_cov_df_long <- strongest_cov_df %>% 
  pivot_longer(cols = c("imbal", "imbal_wt"), names_to = "imbal_type", values_to = "imbal_val") 


p1_full <- p1 + geom_point(data = strongest_cov_df_long[1:num_cov,], 
                           aes(x = imbal_val, y = coeff, z = 0, color = imbal_type), size = psize) + 
  scale_color_manual(labels = c("imbal" = "Pre-wt", "imbal_wt" = "Post-wt"),
                     values = c("imbal" = "red", "imbal_wt" = "forestgreen"))

# teal: #00BFC4
# reddish: #F8766D

p1_full

num_cov_lbl <- 3

hull <- strongest_cov_df[1:(num_cov/2),] %>%
  slice(chull(coeff, imbal))

p1_full_unscaled <- p1_full + 
  ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov_lbl,], 
                            aes(x = imbal, y = coeff, z = 0, label = covar),
                            nudge_y = -0.005, nudge_x = 0.005, label.padding = 0.1,
                            point.padding = 0.1) + 
  # ggrepel::geom_label_repel(data = strongest_cov_df[1:num_cov_lbl,],
  #                           aes(x = imbal_wt, y = coeff, z = 0, label = covar),
  #                           nudge_y = 0.005, nudge_x = 0.005, label.padding = 0.1,
  #                           point.padding = 0.1) + 
  geom_polygon(data = hull, aes(x = imbal, y = coeff, z = 0), alpha = 0.3, fill = "red") + 
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  labs(x = TeX("absolute standardized imbalance in $\\U$"), y = TeX("absolute $\\beta_u$"),
       title = "Disp. Residual (NHANES)", color = c("Imbalance")) + 
  theme(legend.position = "bottom")

# TODO: change font size

p1_full_unscaled



