library(tidyverse)


# load input df

raw_df <- read_csv("../data/data_wide_SM.csv")

ldf <- read_csv("../data/data_SM.csv")


G_var <- c("sex_min" = "SM_inclusive_y_br_ever")

demographics <- c("age" = "age____baseline", "sex" = "sex")


# List 1- based on the ecological system theory
# _ Household: income (or income to needs ratio), 
#.     number of siblings, order of child in family
# _ Family: family conflict/ family history of mental health
# _ Peers: peer victimization
# _ School: school climate measures (e.g., feel safe at school)
# _ Neighborhood: neighborhood safety (perceived), 
#.     geocoded neighborhood measures (area deprivation index)
# _ State: structural stigma against sexual minority

X_list1 <- c("income" = "household_income_mean", "sib_num" = "number_of_siblings____baseline", 
           "sib_order" = "order_among_siblings____baseline", 
           "family_conflict" = "fes_y_ss_fc_mean", "family_mental_health" = "famhx_ss_momdad_scd_p____baseline", 
           "peer_victimization" = "bully_vic_mean", "school_safety" = "srpf_y_ss_ses_mean", 
           "neighborhood_safety" = "neighborhood_crime_y_mean", "adi" = "reshist_addr1_adi_perc_mean", 
           "structural_stigma" = "reshist_state_so_factor_mean")

list1 <- c(G_var, demographics, X_list1)

length(list1)


df1 <- raw_df %>% 
  select(all_of(c(list1))) %>%
  mutate(sib_order = case_when(sib_num == 0 ~ 1,
                               .default = sib_order))


df1 %>% filter(sib_num == 0)


# List 2- based on established risk factors for suicidality in ABCD
# _ Family conflict
# _ Parental monitoring
# _ Weekend screen time
# _ School environment
# _ Negative life events
# _cyberbullying
# _peer victimization
# _racial/ethnic discrimination


