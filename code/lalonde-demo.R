library(Matching)
library(tidyverse)


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


