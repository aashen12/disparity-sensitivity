library(Matching)
library(tidyverse)
library(latex2exp)

andy_ggplot_theme <- function(bs = 20) {
  theme_minimal(base_size = bs) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

data("lalonde")

lalonde <- lalonde %>% 
  mutate(race = ifelse(black == 0 & hisp == 0, "white", "not_white")) %>% 
  tibble()
# mutate with w0, w1, w_rmpw

# A^y: age, married
# N: educ, nodegr

lalonde %>% 
  group_by(race, treat) %>% 
  summarise(mean(re78), mean(educ), mean(nodegr), mean(u74), mean(u75))

lalonde %>% 
  group_by(hisp, treat) %>% 
  summarise(mean(re78), mean(educ), mean(nodegr), mean(u74), mean(u75))

lalonde %>% 
  group_by(black, treat) %>% 
  summarise(mean(re78), mean(educ), mean(nodegr), mean(u74), mean(u75))

lalonde %>% group_by(treat) %>% 
  summarise(mean(re78))

lalonde %>% group_by(race) %>% 
  summarise(mean(re78))

pscore <- glm(black ~ age + married, family = "binomial", data = lalonde)$fitted.values

prob_m_g0_A <- glm(treat ~ age + married, family = "binomial", data = lalonde, weights = 1 - black)$fitted.values

prob_m_g1_AN <- glm(treat ~ age + married + educ + nodegr, family = "binomial", data = lalonde, weights = black)$fitted.values

lalonde <- lalonde %>% mutate(w0 = mean(black == 0) / (1 - pscore), 
                              w1 = mean(black == 1) / pscore, 
                              w_rmpw = w1 * prob_m_g0_A / prob_m_g1_AN)


mu0 <- with(
  lalonde,
  sum(w0 * re78 * (1 - black)) / sum(w0 * (1 - black))
)

mu1 <- with(
  lalonde,
  sum(w1 * re78 * black) / sum(w1 * black)
)

mu_10 <- with(
  lalonde,
  sum(w_rmpw * re78 * black) / sum(w_rmpw * black)
)

mu0
mu1
mu_10

observed <- mu1 - mu0
reduction <- (mu1 - mu_10) |> abs()
residual <- (mu_10 - mu0) |> abs()

observed
reduction
residual

lalonde %>% group_by(black) %>% 
  summarise(me = mean(re78)) %>% 
  mutate(diff = diff(me)) %>% data.frame()

lalonde %>% group_by(treat) %>% 
  summarise(me = mean(re78)) %>% 
  mutate(diff = diff(me)) %>% data.frame()

mu1
mu_10

optBias <- function(R2) {
  corr_bound <- 1 - cor(lalonde$w_rmpw[lalonde$black == 1], lalonde$re78[lalonde$black == 1])^2
  imbalance <- R2 / (1 - R2)
  scaling <- var(lalonde$re78[lalonde$black == 1]) * var(lalonde$w_rmpw[lalonde$black == 1]) 
  return(sqrt(corr_bound * imbalance * scaling))
}

invOptBias <- function(bias) {
  corr_bound <- 1 - cor(lalonde$w_rmpw[lalonde$black == 1], lalonde$re78[lalonde$black == 1])^2
  imbalance <- bias^2 / (corr_bound * var(lalonde$re78[lalonde$black == 1]) * var(lalonde$w_rmpw[lalonde$black == 1]))
  return(imbalance / (1 + imbalance))
}


R2_reduction <- invOptBias(reduction) 
print(
  paste0("If an unmeasured confounder explains at least ", 
         round(R2_reduction * 100, 2), "% of the variation in the true weights, then the disparity reduction becomes zero.")
) 

R2_residual <- invOptBias(residual)
print(
  paste0("If an unmeasured confounder explains at least ", 
         round(R2_residual * 100, 2), "% of the variation in the true weights, then the disparity remaining becomes zero.")
) 

R2 <- c(seq(0.1, 0.6, by = 0.1), R2_residual, 0.7)
bias <- optBias(R2)
# lower_pe_reduc <- reduction - bias
# upper_pe_reduc <- reduction + bias
lower_pe_resid <- residual - bias
upper_pe_resid <- residual + bias

df_resid <- data.frame(R2 = round(R2, 3), bias, lower_pe_resid, upper_pe_resid) 


# make plot with using ggplot with R2 on x-axis (discretized) and estimated residual on y-axis with error bars for lower and upper bounds
df_resid %>% 
  ggplot(aes(x = as.character(R2))) + 
  geom_errorbar(aes(ymin = lower_pe_resid, ymax = upper_pe_resid), color = "dodgerblue2", 
                width = 0.5, linewidth = 2, alpha = 0.8) + 
  geom_point(y = residual, color = "black", size = 6) +
  geom_hline(yintercept = 0, linewidth = 1, alpha = 0.5) + 
  labs(x = expression(R^2), y = TeX(r'($\hat{\mu}_1^0 - \hat{\mu}_0$)'), title = "Residual Disparity") + 
  andy_ggplot_theme(22)
  




R2 <- seq(R2_reduction, 0.1, by = 0.02)
bias <- optBias(R2)
lower_pe_reduc <- reduction - bias
upper_pe_reduc <- reduction + bias


df_reduc <- data.frame(R2 = round(R2, 3), bias, lower_pe_reduc, upper_pe_reduc)

# make plot with using ggplot with R2 on x-axis (discretized) and estimated reduction on y-axis with error bars for lower and upper bounds

df_reduc %>% 
  ggplot(aes(x = as.character(R2))) + 
  geom_errorbar(aes(ymin = lower_pe_reduc, ymax = upper_pe_reduc), color = "firebrick1", 
                width = 0.5, linewidth = 1.7, alpha = 0.8) + 
  geom_point(y = reduction, color = "black", size = 3.5) +
  geom_hline(yintercept = 0, linewidth = 1, alpha = 0.5) + 
  labs(x = TeX(r'($R^2$)'), y = TeX(r'($\hat{\mu}_1 - \hat{\mu}_1^0$)'), title = "Disparity Reduction") + 
  andy_ggplot_theme(22)
  





