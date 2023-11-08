library(Matching)
library(tidyverse)

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

pscore <- glm(black ~ age + married, family = "binomial", data = lalonde)$fitted.values

prob_m_g0_A <- glm(treat ~ age + married, family = "binomial", data = lalonde, weights = 1 - black)$fitted.values

prob_m_g1_AN <- glm(treat ~ age + married + educ + nodegr, family = "binomial", data = lalonde, weights = black)$fitted.values

lalonde <- lalonde %>% mutate(w0 = mean(black == 0) / (1 - pscore), 
                              w1 = mean(black == 1) / pscore, 
                              w_rmpw = w1 * prob_m_g1_AN / prob_m_g0_A)

View(lalonde)

mu_10 <- 

