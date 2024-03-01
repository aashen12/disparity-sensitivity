# perform analysis of add health data
# devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))

library(tidyverse)
library(decompsens)
library(latex2exp)

set.seed(122)

Hajek <- function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}

HT <- function(x, w, na.rm = TRUE) {
  mean(x * w, na.rm = na.rm)
}


df_yz <- read_csv("../../data/add-health-data/trt_out_df.csv") %>% select(-attempt) %>% drop_na()
df_x <- read_csv("../../data/add-health-data/covariate_df.csv") %>% select(-dad_ed) %>% drop_na()

covariates <- names(df_x)[!names(df_x) %in% c("AID", "CLUSTER2", "GSWGT1")]

allowables <- c("age", "sex")

non_allowables <- setdiff(covariates, allowables)

df_full <- inner_join(df_yz, df_x, by = "AID")
# G_tot <- df_full$sex_minority
# G1_ind <- which(G_tot == 1)
# G0_ind <- sample(which(G_tot == 0), size = sum(G_tot == 1), replace = TRUE)
# 
# df <- df_full %>% slice(c(G1_ind, G0_ind))

df <- df_full

XA <- model.matrix(~ . -1, data = df %>% select(all_of(allowables)))
XN <- model.matrix(~ . -1, data = df %>% select(all_of(non_allowables)))

X <- cbind(XA, XN)

head(X)

Y <- df$ideation
Z <- df$school_belong
G <- df$sex_minority

mu1 <- mean(Y[G == 1])
mu0 <- mean(Y[G == 0])



obs_disp <- mu1 - mu0
obs_disp

sd_obs_disp <- sqrt(var(Y[G == 1]) / sum(G == 1) + var(Y[G == 0]) / sum(G == 0))

# construct 95% CI for obs_disp
ci <- c(obs_disp - qnorm(0.95) * sd_obs_disp, obs_disp + qnorm(0.95) * sd_obs_disp)
ci


mu1/mu0

trim <- 0
allow <- F

w <- decompsens::estimateRMPW(G=G, Z=Z, Y=Y, XA=XA, XN=XN, trim = trim, allowable = allow)
summary(w)

mu10 <- sum(Y[G == 1] * w[G == 1]) / sum(w[G == 1])
mu10
mu1
mu0

mu1-mu0


reduction <- mu1 - mu10
reduction
mu1/mu10

residual <- mu10 - mu0
residual
mu10/mu0

## Bootstrap standard errors
B <- 2000
resid_boot <- numeric(B)
out <- parallel::mclapply(1:B, function(i) {
  ind <- sample(1:length(Y), length(Y), replace = TRUE)
  w_boot <- decompsens::estimateRMPW(G=G[ind], Z=Z[ind], Y=Y[ind], XA=XA[ind,], XN=XN[ind,], trim = trim, allowable = allow)
  mu10_boot <- sum(Y[ind][G[ind] == 1] * w_boot[G[ind] == 1]) / sum(w_boot[G[ind] == 1])
  mu1_boot <- mean(Y[ind][G[ind] == 1])
  mu0_boot <- mean(Y[ind][G[ind] == 0])
  resid_boot <- mu10_boot - mu0_boot
  red_boot <- mu1_boot - mu10_boot
  list(resid_boot = resid_boot, red_boot = red_boot)
})

out_resid <- unlist(lapply(out, function(x) x[["resid_boot"]]))
out_red <- unlist(lapply(out, function(x) x[["red_boot"]]))

se_resid <- sqrt((1 / (B-1)) * sum((out_resid - mean(out_resid))^2))
se_resid
# 90% CI
ci_resid <- c(residual - qnorm(0.95) * se_resid, residual + qnorm(0.95) * se_resid)
ci_resid

se_red <- sqrt((1 / (B-1)) * sum((out_red - mean(out_red))^2))
se_red
# 90% CI
ci_red <- c(reduction - qnorm(0.8) * se_red, reduction + qnorm(0.8) * se_red)
ci_red


###########################################################################################
###########################################################################################
###########################################################################################

### Disparity Residual ###

estimand <- "res"

extrema <- decompsens::getExtrema(G=G, Y=Y, w=w, gamma = log(1), estimand = estimand, RD = T)
extrema

boot_ci <- decompsens::bootstrapCI(G=G, Z=Z, Y=Y, XA=XA, XN=XN, w=w, gamma = log(1.04),
                                   alpha=0.1, B=2000, estimand = estimand, stratify = FALSE,
                                   parallel = TRUE, allowable = allow, RD = T)
boot_ci

lam_resid <- 1.05

# cis <- sapply(
#   Lam_grid,
#   function(lam) {
#     decompsens::bootstrapCI(G=G, Z=Z, Y=Y, XA=XA, XN=XN, w=w, gamma = lam,
#                             alpha=0.05, B=2000, estimand = estimand)
#   }
# )
# 
# if (ncol(cis) == length(Lam_grid)) {
#   nlam <- length(Lam_grid)
# } else {
#   print("Error in bootstrapping")
# }
# 
# cis <- cis %>% t()

Lam <- lam_resid

bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, Lambda = Lam)
amplification <- decompsens::informalAmplify(G, Z, XA, XN, Y, Lambda = Lam)

strongest_cov_df <- amplification[[1]] %>% drop_na()
max_imbal_stnd <- amplification$max_imbal_stnd
max_betau_01 <- amplification$max_beta
str(amplification)

maxbias <- max(abs(bounds[[1]]))

max_betau_01 <- max(strongest_cov_df$coeff)
max_imbal <- max(strongest_cov_df$imbal)

beta <- seq(min(strongest_cov_df$coeff) - 0.02, max(strongest_cov_df$coeff) + 0.02, by = 0.01)
imbalance <- seq(0.001, max(strongest_cov_df$imbal) + 0.2, by = 0.005)
data.fit <- expand.grid(beta, imbalance)
names(data.fit) <- c("beta", "imbalance")

df_plot <- data.fit %>% 
  mutate(bias = abs(beta * imbalance))

summary(df_plot$bias)

round_down <- function(a) {
  floor(a * 100) / 100
}


bins <- seq(max(df_plot$bias) - 1.6 * sd(df_plot$bias),  max(df_plot$bias) + 1 * sd(df_plot$bias), length.out = 4) %>% 
  round_down()

num_cov <- 10 # specify 2x what you want
psize <- 6.5

p1 <- df_plot %>%
  ggplot(aes(x = imbalance, y = beta, z = bias)) +
  theme_bw() + 
  geom_contour(col="gray55", breaks = bins) + 
  metR::geom_text_contour(aes(z = bias), 
                          breaks = bins,
                          stroke = 0.2, skip = 0) + 
  metR::geom_contour_fill(breaks = c(maxbias, 1000 * maxbias), fill='powderblue', alpha = 0.5) +
  geom_contour(breaks = c(maxbias), col='blue', linewidth = 1) + 
  metR::geom_text_contour(aes(z = maxbias), stroke = 0.2) + 
  geom_point(x = 0.2, y = maxbias / 0.2, size = psize - 0.2, color = "black") +
  geom_text(x = 0.2 + 0.02 , y = maxbias / 0.2 + 0.004, 
            label = paste0("Max Bias: ", round(maxbias, 3)), size = psize, color = "black")
  
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


with(strongest_cov_df, chull(coeff, imbal))

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
       title = "Disp. Residual (Add Health)", color = c("Imbalance")) + 
  theme(legend.position = "bottom")

# TODO: change font size

p1_full_unscaled
