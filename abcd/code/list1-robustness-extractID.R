# Doing a "robustness check" for temporality
# only considers kids who do no t have suicide ideation at baseline
rm(list = ls())
library(tidyverse)

set.seed(122357)
# load input df

raw_df <- read_csv("../data/data_SM.csv")
head(raw_df)

unique(raw_df$eventname)

df_baseline <- raw_df %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  select(id = src_subject_id, si = SI_y) %>% 
  filter(si == 1)

si_baseline_id <- df_baseline$id

saveRDS(si_baseline_id, "../data/si_baseline_ids.rds")
