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

lalonde %>% group_by(treat, black) %>% 
  summarise(mean(re78))

lalonde %>% group_by(race) %>% 
  summarise(mean(re78))


e0 <- glm(treat ~ age + married, family = "binomial", data = lalonde, weights = 1 - black)$fitted.values
e1 <- glm(treat ~ age + married + educ + nodegr, family = "binomial", data = lalonde, weights = black)$fitted.values

lalonde <- lalonde %>% mutate(
  w1 = e0 / e1,
  w0 = (1 - e0) / (1 - e1),
  w_rmpw = ifelse(treat == 1, w1, w0)
)

mu0 <- with(
  lalonde,
  sum(re78 * (1 - black)) / sum(1 - black)
)

mu1 <- with(
  lalonde,
  sum(re78 * black) / sum(black)
)

mu_10 <- with(
  lalonde,
  sum(w_rmpw * re78 * black) / sum(w_rmpw * black)
)

mu0
mu1
mu_10


cf1 <- with(
  lalonde,
  sum(w_rmpw * re78 * black) / sum(w_rmpw * black)
)

cf2 <- with(
  lalonde, 
  {
    num <- sum(black * treat * re78 * w1) + sum(black * (1 - treat) * re78 * w0)
    den <- sum(black * treat * w1) + sum(black * (1 - treat)  * w0)
    num / den
  }
)

cf3 <- with(
  lalonde, 
  {
    treat <- sum(black * treat * re78 * w1) / sum(black * treat * w1)
    control <-  sum(black * (1 - treat) * re78 * w0) / sum(black * (1 - treat) * w0)
    treat + control
  }
)

cf1
cf2
cf3


## Total disparity

tot_disp <- lm(re78 ~ black, data = lalonde)$coef[2]
mu1 - mu0


##################################################################
##                      Residual Disparity                      ##
##################################################################

disp_residual <- lm(re78 ~ black, data = lalonde, weights = ifelse(black == 1, w_rmpw, 1))$coef[2]
disp_residual
c(cf1, cf2, cf3) - mu0

#################################################################
##                     Disparity Reduction                     ##
#################################################################

df_black <- lalonde %>% filter(black == 1)

df_reduction <- bind_rows(
  df_black %>% mutate(D = 1),
  df_black %>% mutate(D = 0)
)

disp_reduction <- lm(re78 ~ D, data = df_reduction, weights = ifelse(D == 0, w_rmpw, 1))$coef[2]
disp_reduction
mu1 - c(cf1, cf2, cf3)



# checking observed disparity values mu1 and mu0
lalonde %>% group_by(black) %>% 
  summarise(me = mean(re78)) %>% 
  mutate(diff = diff(me)) %>% data.frame()

mu1
mu0

# test
