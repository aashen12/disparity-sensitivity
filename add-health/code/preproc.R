# Pre-processing Add Health dataset

library(tidyverse)
library(prettyR)
library(bannerCommenter)

coerceResponse <- function(MYVAR) {
  # Function to remove unnecessary labels from response entries
  # essentially 'cleans' the entries
  lbls <- sort(levels(MYVAR))
  lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
  out <- add.value.labels(
    as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", MYVAR)), 
    lbls
  )
  out
}



#################################################################
##         Loading data and identifying question codes         ##
#################################################################

load("../../data/add-health-data/DS1.rda") # respondent answer
load("../../data/add-health-data/DS2.rda") # parent answer
load("../../data/add-health-data/DS3.rda")
load("../../data/add-health-data/DS4.rda")

ds1 <- da21600.0001
ds2 <- da21600.0002
ds3 <- da21600.0003
ds4 <- da21600.0004

demographic_codes <- c(
  latinx = "H1GI4", # Hisp/Latino
  white = "H1GI6A", # White
  black = "H1GI6B", # Black
  native = "H1GI6C", # Native American
  asian = "H1GI6D", # Asian
  sex = "BIO_SEX" # Gender, 1 = male, 2 = female
)

# alcohol_code <- "H1TO" # section 28
# sex_drunk_code <- "H1JO" 
# contraception_code <- "H1CO" # section 24
# knowledge_quiz_code <- "H1KQ1A"
# alc_in_home_code <- "H1TO51"

treatment_codes <- c(
  parent_care = "H1PR3",
  adult_care = "H1PR1",
  psych_counsel = "H1HS3", # In the past year, have you received psychological or emotional counseling?
  phys_exam = "H1HS1", # In the past year have you had a routine physical examination?
  med_care = "H1GH26", # you thought you should get medical care, but you did not?
  suicide_educ = "H1TS17", # received suicide educ
  school_belong = "H1ED20", # belonging in school
  school_friends = "H1ED19" # friends in school
)

sexuality_codes <- c(
  attract_fem = "H1NR1", # attracted to female?
  attract_male = "H1NR2" # attracted to male
)


outcome_codes <- c(
  ideation = "H1SU1",
  attempt = "H1SU2"
)



smoke_household <- "PA63"


codes_y_z <- c("AID", demographic_codes, treatment_codes, outcome_codes, sexuality_codes)


df_raw_y_z <- ds1 %>% 
  dplyr::select(all_of(codes_y_z)) %>% 
  inner_join(ds4, by = "AID")


df_yz <- apply(df_raw_y_z, 2, coerceResponse) %>% data.frame() %>% tibble() 

table(df_yz$school_belong)
table(df_yz$school_friends)


# Assuming 'df' is your data frame with dummy variables
# and the dummy variable columns are 'latinex', 'white', 'black', 'native', 'asian'
race_columns <- c('latinx', 'white', 'black', 'native', 'asian')

# Find the index of the max value in each row (which dummy variable is 1)
race_index <- max.col(df_yz[race_columns], ties.method = "first")

# Map the index to the actual race names
df_yz$race <- race_columns[race_index] #%>% factor()

df_yz <- df_yz %>% 
  mutate_at(vars("attempt"), 
            ~case_when(.x == 0 ~ 0, .x > 0 ~ 1, .default = .x)) %>%
  mutate_at(vars("school_belong", "school_friends"), 
            ~case_when(
              .x %in% c(1, 2) ~ 1,
              .x > 2 ~ 0,
            )) %>%
  mutate_at(vars("AID"), factor) %>% 
  mutate_at(vars("sex"), ~case_when(.x == 1 ~ "M", .x == 2 ~ "F")) %>% 
  mutate_at(vars("sex"), factor) %>% 
  mutate(
    sex_minority = case_when(
      sex == "F" & attract_fem == 1 ~ 1,
      sex == "M" & attract_male == 1 ~ 1,
      .default = 0
    ) %>% factor()
  ) %>% select(-attract_fem, -attract_male)

# Now df$race contains the categorical variable for race

races <- c("latinx", "asian")



Hajek <- function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}

HT <- function(x, w, na.rm = TRUE) {
  mean(x * w, na.rm = na.rm)
}

df_yz %>% 
  group_by(sex_minority) %>%
  summarise(mean_attempt = Hajek(attempt, w = GSWGT1, na.rm = TRUE),
            mean_ideation = Hajek(ideation, w = GSWGT1, na.rm = TRUE),
            sd_ideation = sd(ideation, na.rm = TRUE),
            n = n())

df_yz %>% 
  group_by(sex_minority) %>%
  summarise(mean_adult = Hajek(adult_care, w = GSWGT1, na.rm = TRUE),
            mean_parent = Hajek(parent_care, w = GSWGT1, na.rm = TRUE),
            n = n())


df_yz %>% 
  group_by(sex_minority) %>%
  summarise(mean_psych = Hajek(psych_counsel, w = GSWGT1, na.rm = TRUE),
            mean_phys = Hajek(phys_exam, w = GSWGT1, na.rm = TRUE),
            mean_med = Hajek(med_care, w = GSWGT1, na.rm = TRUE),
            n = n())

df_yz %>%
  group_by(sex_minority) %>%
  summarise(mean_suicide_educ = Hajek(suicide_educ, w = GSWGT1, na.rm = TRUE))

df_yz %>% 
  group_by(sex_minority) %>% 
  summarise(mean_school_belong = Hajek(school_belong, w = GSWGT1, na.rm = TRUE),
            mean_school_friends = Hajek(school_friends, w = GSWGT1, na.rm = TRUE))


trt_out_df <- df_yz %>% 
  dplyr::select(AID, sex_minority, school_belong, attempt, ideation, weight = GSWGT1)

write_csv(trt_out_df, file = "../data/trt_out_df.csv")
# outcome: suicide
# treatment: belonging in school

## Now we do the covariates ##

## Baseline: age, sex, race, household income, SES, parental education/marital status

# FROM RAN: List 1- based on the ecological system theory, where we choose one-two exposures from each level of exposure:
#   _ Household: income (or income to needs ratio)
# 
# _ Family: family conflict/ family history of mental health
# 
# _ School: school climate measures
# 
# _ Neighborhood: neighborhood safety (perceived), geocoded neighborhood measures (area deprivation index)
# 
# _ State: structural stigma against sexual minority

# List 2- based on established risk factors for suicidality in ABCD (based on others' and our past works)
# 
# _ Family conflict
# 
# _ Parental monitoring
# 
# _ Weekend screen time
# 
# _ School environment
# 
# _ Negative life events

baseline_covs <- c(
  sex = "BIO_SEX", # Gender, 1 = male, 2 = female
  income = "PA55", # numeric number 
  mom_ed = "H1RM1", # number for increasing level
  dad_ed = "H1RF1", # number for increasing level
  school_safe = "H1ED24" # do u feel safe at school 1: agree - 5: disagree
)

age = "S1" # separate since we do not want to select contains("S1")

list1_covs <- c(
  neigh_safety = "H1NB5" # do u feel safe in neigh? 0 no 1 yes
)

# list2_covs <- c(
#   
# )

codes_x <- c("AID", baseline_covs, list1_covs)

race_df <- df_yz %>% select(AID, race) %>% mutate_at("race", factor)

df_raw_x <- ds1 %>% 
  data.frame() %>% tibble() %>% 
  mutate_at("AID", factor) %>% 
  dplyr::select(all_of(codes_x), age = S1) %>% 
  inner_join(ds4, by = "AID") 

df_x <- apply(df_raw_x, 2, coerceResponse) %>% data.frame() %>% tibble() %>% 
  mutate_at("AID", factor) %>% 
  inner_join(race_df, by = "AID") %>% 
  relocate(race, .after = "AID") %>% 
  mutate_at("sex", ~case_when(.x == 1 ~ "M", .x == 2 ~ "F") %>% factor()) %>% 
  relocate(age, .after = "AID") %>% 
  mutate_at(vars("mom_ed", "dad_ed"), 
            ~case_when(
              .x == 10 ~ 1,
              .x == 11 ~ NA,
              .x == 12 ~ NA,
              .default = .x))

write_csv(df_x, file = "../data/covariate_df.csv")


