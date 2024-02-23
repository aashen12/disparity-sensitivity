##############
# Choose data
##############
# options: 'rhc', 'nhanes', 'lalonde', 'school', 'fish', 'lalonde interactions'
data <- 'fish'

# save inputs into input data frame 
input_df <-
  data.frame(
    data = data,
    stringsAsFactors = F
  )

# save
saveRDS(input_df, file = "../inputs/input_df.rds")

