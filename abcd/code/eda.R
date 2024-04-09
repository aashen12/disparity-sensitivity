library(tidyverse)

set.seed(122357)
# load input df

raw_df <- read_csv("../data/data_wide_SM.csv")

ldf <- read_csv("../data/data_SM.csv")

## Construct binarized version of parental acceptance ##

## Parental acceptance: crpbi_y_ss_parent_mean crpbi_y_ss_caregiver_mean

parent <- raw_df$crpbi_y_ss_parent_mean
caregiver <- raw_df$crpbi_y_ss_caregiver_mean
sm <- raw_df$SM_inclusive_y_br_ever

summary(parent)
summary(caregiver)


data.frame(parent, caregiver, sm) %>% 
  filter(!is.na(parent) & !is.na(caregiver)) %>% 
  pivot_longer(cols = c(parent, caregiver), names_to = "role", values_to = "score") %>%
  ggplot(aes(x = score, fill = role)) +
  geom_density(color = "black", alpha = 0.5) +
  labs(x = "score", y = "density", title = "Density plot of parental and caregiver acceptance scores",
       fill = "Role") + 
  facet_wrap(~sm) +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  theme(text = element_text(size = 20)) 



accept <- mean(parent + caregiver, na.rm = TRUE)

hist(parent + caregiver)

# Criteria for 1:
#. Scores of 3 for both
#. Scores of 3 for parent (if no caregiver)
#. Scores of 3 for both OR score of 3 for parent and at least 2 for caregiver

good_parent_thresh <- 8/3

df_withz <- raw_df %>% 
  mutate(parent_accept = case_when(
    #crpbi_y_ss_parent_mean >= good_parent_thresh & crpbi_y_ss_caregiver_mean >= good_parent_thresh ~ 1,
    crpbi_y_ss_parent_mean == 3 & crpbi_y_ss_caregiver_mean == 3 ~ 1,
    crpbi_y_ss_parent_mean == 3 & is.na(crpbi_y_ss_caregiver_mean) ~ 1,
    crpbi_y_ss_parent_mean == 3 & crpbi_y_ss_caregiver_mean > good_parent_thresh ~ 1,
    crpbi_y_ss_parent_mean > good_parent_thresh & crpbi_y_ss_caregiver_mean == 3 ~ 1,
    is.na(crpbi_y_ss_parent_mean) & crpbi_y_ss_caregiver_mean == 3 ~ 1,
    crpbi_y_ss_parent_mean < 3 | crpbi_y_ss_caregiver_mean < 3 ~ 0,
    .default = NA
  ))


df_withz$parent_accept %>% table()


## Basics to adjust for and include in the model ##


Z_var <- c("parent_accept" = "parent_accept")


Y_var <- c("ideation" = "SI_y_ever") # ideation
#Y_var <- c("ideation" = "SA_y_ever") # attempt

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
           "peer_victimization" = "bully_vic_mean", 
           "school_safety" = "srpf_y_ss_ses_mean", 
           "neighborhood_safety" = "neighborhood_crime_y_mean", "adi" = "reshist_addr1_adi_perc_mean", 
           "structural_stigma" = "reshist_state_so_factor_mean")

list1 <- c(G_var, demographics, X_list1, Y_var, Z_var)

length(list1)


df1 <- df_withz %>% 
  select(all_of(list1)) %>%
  mutate(sib_order = case_when(sib_num == 0 ~ 1,
                               .default = sib_order)) %>% 
  drop_na(sex_min, ideation, parent_accept) #%>% drop_na()


df1 %>% group_by(sex_min) %>% 
  summarise(n = n(), mean(ideation), mean(parent_accept))


df_x1 <- df1 %>% select(all_of(names(c(demographics, X_list1))))
write_csv(df_x1, "../data/list1_X.csv")


df_yz1 <- df1 %>% select(all_of(c(names(Y_var), names(Z_var), names(G_var), "peer_victimization")))
write_csv(df_yz1, "../data/list1_YZG.csv")

Y <- df1$ideation
G <- df1$sex_min
Z <- df1$parent_accept

table(Z)
table(Y)
table(G)
any(is.na(Z))



# List 2- based on established risk factors for suicidality in ABCD
# _ Family conflict
# _ Parental monitoring
# _ Weekend screen time
# _ School environment
# _ Negative life events
# _cyberbullying
# _peer victimization
# _racial/ethnic discrimination

X_list2 <- c("family_conflict" = "fes_y_ss_fc_mean", 
             "parental_monitoring" = "fes_y_ss_pm_mean", 
             "weekend_screen_time" = "srpf_y_ss_wst_mean", 
             "school_environment" = "srpf_y_ss_ses_mean", 
             "negative_life_events" = "fes_y_ss_nle_mean", 
             "cyberbullying" = "bully_cyb_mean", 
             "peer_victimization" = "bully_vic_mean")






