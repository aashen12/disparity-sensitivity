library(Matching)
library(tidyverse)

data("lalonde")

lalonde <- lalonde %>% 
  mutate(race = ifelse(black == 0 & hisp == 0, "white", "not_white")) 
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




