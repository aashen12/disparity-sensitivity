rm(list = ls())
library(tidyverse)

set.seed(122357)
# load input df

raw_df <- read_csv("../data/data_wide_SM.csv")

# ldf <- read_csv("../data/data_SM.csv")

Z_method <- 6 # aggregate = 0, or specify number 1-5. if 6 then we aggregate the better_worry + better_upset
Z_method_questions <- c(
  "better_worry",
  "smile",
  "better_upset",
  "love",
  "easy_talk"
)

##################################################################
##           Parental and caregiver acceptance scores           ##
##################################################################

# crpbi_parent1_y	Makes me feel better after talking over my worries with him/her (1)
# crpbi_parent2_y Smiles at me very often. (2)
# crpbi_parent3_y Is able to make me feel better when I am upset. (3)
# crpbi_parent4_y Believes in showing his/her love for me. (4)
# crpbi_parent5_y Is easy to talk to. (5)
# crpbi_caregiver1_y Is there a second adult who cares for you, that you spend a significant amount of time with, like your other parent, step-parent, grandparent, aunt or uncle?
# crpbi_caregiver2_y	Who is the second caregiver? If more than one caregiver, choose the person you spend the most time with.
# crpbi_caregiver12_y	Makes me feel better after talking over my worries with them.
# crpbi_caregiver13_y	Smiles at me very often.
# crpbi_caregiver14_y	Is able to make me feel better when I am upset.
# crpbi_caregiver15_y	Believes in showing their love for me.
# crpbi_caregiver16_y	Is easy to talk to.




if (Z_method == 0) {
  thresh <- 0
  Z_method_question <- "aggregate"
  parent_questions <- paste0("crpbi_parent", 1:5, "_y_mean")
  caregiver_questions <- paste0("crpbi_caregiver", 12:16, "_y_mean") # crpbi_caregiver12_y
  df_withz <- raw_df %>% 
    #select(src_subject_id, all_of(parent_questions), all_of(caregiver_questions)) %>% 
    mutate(across(all_of(parent_questions), 
                  ~ case_when(. >= 2.99 ~ 1, 
                              . < 2.99 ~ 0,
                              .default = NA))) %>% 
    mutate(across(all_of(caregiver_questions),
                  ~ case_when(. >= 2.99 ~ 1, 
                              . < 2.99 ~ 0,
                              .default = NA))) %>% # binarize: if 3, set to 1. If < 3, set to 0. This "3" may be replaced with 3 - \epsilon
    mutate(prop_1_parent = rowMeans(select(., all_of(parent_questions)), na.rm = FALSE),
           prop_1_caregiver = rowMeans(select(., all_of(caregiver_questions)), na.rm = FALSE), 
           .before = crpbi_parent1_y_mean) %>% # compute proportion of 3s for each row
    mutate(parent_accept = case_when(
      prop_1_parent == 1 & prop_1_caregiver >= 1 - thresh ~ 1,
      prop_1_parent >= 1 - thresh & prop_1_caregiver == 1 ~ 1,
      prop_1_parent == 1 & is.na(prop_1_caregiver) ~ 1,
      prop_1_parent == 0 & prop_1_caregiver <= thresh ~ 0,
      prop_1_parent <= thresh & prop_1_caregiver == 0 ~ 0,
      prop_1_parent == 0 & is.na(prop_1_caregiver) ~ 0
    ))
  print(paste0("Z_method = ", Z_method))
} else if (Z_method %in% 1:5) {
  Z_method_question <- Z_method_questions[Z_method]
  parent_questions <- paste0("crpbi_parent", 1:5, "_y_mean")
  caregiver_questions <- paste0("crpbi_caregiver", 12:16, "_y_mean") # crpbi_caregiver12_y
  q_index <- Z_method
  parent_column <- paste0("crpbi_parent", q_index, "_y_mean")
  caregiver_column <- paste0("crpbi_caregiver", q_index + 11, "_y_mean")
  df_withz <- raw_df %>%
    rename(prop_1_parent = all_of(parent_column), prop_1_caregiver = all_of(caregiver_column)) %>%
    #select(src_subject_id, prop_1_parent = all_of(parent_column), prop_1_caregiver = all_of(caregiver_column)) %>%
    mutate(prop_1_parent = ifelse(is.na(prop_1_parent), NA, 
                                  ifelse(prop_1_parent == max(prop_1_parent, na.rm = TRUE), 1, 0))) %>%
    mutate(prop_1_caregiver = ifelse(is.na(prop_1_caregiver), NA, 
                                     ifelse(prop_1_caregiver == max(prop_1_caregiver, na.rm = TRUE), 1, 0))) %>%
    mutate(parent_accept = case_when(
      prop_1_parent == 1 & prop_1_caregiver == 1 ~ 1,
      prop_1_parent == 1 & is.na(prop_1_caregiver) ~ 1,
      prop_1_parent == 0 & prop_1_caregiver == 0 ~ 0,
      prop_1_parent == 0 & is.na(prop_1_caregiver) ~ 0
    ))
  print(paste0("Z_method = ", Z_method))
} else if (Z_method == 6) {
  thresh <- 3 + 3 + 3 + 2.99#2.65
  Z_method_question <- "worry_upset"
  parent_column <- paste0("crpbi_parent", c(1, 3), "_y_mean")
  caregiver_column <- paste0("crpbi_caregiver", c(1, 3) + 11, "_y_mean")
  df_withz <- raw_df %>% 
    rename(parent_worry = all_of(parent_column[1]), parent_upset = all_of(parent_column[2]),
           caregiver_worry = all_of(caregiver_column[1]), caregiver_upset = all_of(caregiver_column[2])) %>%
    #select(src_subject_id, parent_worry, parent_upset, caregiver_worry, caregiver_upset) %>% 
    group_by(src_subject_id) %>%
    mutate(parent_accept = case_when(
      sum(parent_worry, parent_upset, caregiver_worry, caregiver_upset, na.rm = FALSE) > thresh ~ 1,
      !is.na(caregiver_worry) & !is.na(caregiver_upset) & all(parent_worry <= 2.99, parent_upset <= 2.99, caregiver_worry <= 2.99, caregiver_upset <= 2.99) ~ 0,
      sum(parent_worry, parent_upset) >= 5.99 & is.na(caregiver_worry) & is.na(caregiver_upset) ~ 1,
      all(parent_worry < 2.99, parent_upset < 2.99) & is.na(caregiver_worry) & is.na(caregiver_upset) ~ 0,
      .default = NA
    )) %>% ungroup()
    print(paste0("Z_method = ", Z_method))
}

Z_var <- c("parent_accept" = "parent_accept")

# NDAR_INVXMJE5DN0




Y_var <- c("ideation" = "SI_y_ever") # ideation
#Y_var <- c("attempt" = "SA_y_ever") # attempt

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
  select(src_subject_id, all_of(list1)) %>%
  mutate(sib_order = case_when(sib_num == 0 ~ 1,
                               .default = sib_order)) %>% 
  drop_na(sex_min, ideation, parent_accept) #%>% drop_na()


df1 %>% group_by(sex_min) %>% 
  summarise(n = n(), mean(ideation), mean(parent_accept))


df_x1 <- df1 %>% select(src_subject_id, all_of(names(c(demographics, X_list1))))
write_csv(df_x1, paste0("../data/list1_X_", Z_method_question, ".csv"))


df_yz1 <- df1 %>% select(src_subject_id, all_of(c(names(Y_var), names(Z_var), names(G_var), "peer_victimization")))
write_csv(df_yz1, paste0("../data/list1_YZG_", Z_method_question, ".csv"))

message(paste0("Wrote df for Z method: ", Z_method_question))

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

##################################################################
##                           OLD CODE                           ##
##################################################################

# Criteria for 1:
#. Scores of 3 for both
#. Scores of 3 for parent (if no caregiver)
#. Scores of 3 for both OR score of 3 for parent and at least 2 for caregiver

# good_parent_thresh <- 8/3
# 
# df_withz <- raw_df %>% 
#   mutate(parent_accept = case_when(
#     #crpbi_y_ss_parent_mean >= good_parent_thresh & crpbi_y_ss_caregiver_mean >= good_parent_thresh ~ 1,
#     crpbi_y_ss_parent_mean == 3 & crpbi_y_ss_caregiver_mean == 3 ~ 1,
#     crpbi_y_ss_parent_mean == 3 & is.na(crpbi_y_ss_caregiver_mean) ~ 1,
#     crpbi_y_ss_parent_mean == 3 & crpbi_y_ss_caregiver_mean > good_parent_thresh ~ 1,
#     crpbi_y_ss_parent_mean > good_parent_thresh & crpbi_y_ss_caregiver_mean == 3 ~ 1,
#     is.na(crpbi_y_ss_parent_mean) & crpbi_y_ss_caregiver_mean == 3 ~ 1,
#     crpbi_y_ss_parent_mean < 3 | crpbi_y_ss_caregiver_mean < 3 ~ 0,
#     .default = NA
#   ))


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

