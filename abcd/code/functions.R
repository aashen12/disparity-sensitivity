NAImpute <- function(df, append = FALSE) {
  for(i in c(1:ncol(df))){
    if(any(is.na(df[,i]))){
      print(paste('Missing values found in column', i ,'of X; imputing and adding missingness indicators'))
      if (append) {
        df <- cbind(df, is.na(df[,i]))
        colnames(df)[ncol(df)] <- paste(colnames(df)[i],'NA', sep = '')
      }
      if(all(df[,i][!is.na(df[,i])] %in% c(0, 1))) {
        df[which(is.na(df[,i])),i] <- names(sort(-table(df[,i])))[1] %>% as.numeric()
        df[which(is.na(df[,i])),i] <- as.numeric(df[which(is.na(df[,i])),i])
      } else {
        df[which(is.na(df[,i])),i] <- mean(df[,i], na.rm = TRUE)
      }
    }
  }
  df
}



create_model_matrix <- function(df) {
  library(splines)
  
  # Identify continuous and binary variables
  binary_vars <- names(df)[sapply(df, function(x) all(x %in% c(0, 1)) && is.numeric(x))]
  continuous_vars <- setdiff(names(df), binary_vars)
  
  # Start formula string with original variables
  formula_terms <- c(continuous_vars, binary_vars)
  
  # Add natural cubic splines for continuous variables
  formula_terms <- c(formula_terms, 
                     paste0("ns(", continuous_vars, ", df = 3)"))
  
  # Add 2-way interactions for binary variables
  if (length(binary_vars) > 1) {
    interaction_terms <- apply(combn(binary_vars, 2), 2, 
                               function(x) paste0(x[1], ":", x[2]))
  } else {
    interaction_terms <- NULL
  }
  formula_str <- paste(c(formula_terms, interaction_terms), collapse = " + ")
  
  # Create the model matrix
  model_mat <- model.matrix(as.formula(paste0("~ 0 + ", formula_str)), data = df)
  
  return(model_mat)
}



lovePlot <- function(pre_weight, post_weight, num_covs = 100, title = "Covariate Balance") {
  subset_covs <- 2 * min(num_covs, length(pre_weight))
  pre_weight <- abs(pre_weight)
  post_weight <- abs(post_weight)
  if (any(names(pre_weight) != names(post_weight))) {
    stop("Pre and post weights must have the same names")
  }
  data.frame(
    covariate = names(pre_weight),
    pre_weight = pre_weight,
    post_weight = post_weight
  ) %>% 
    arrange(desc(pre_weight)) %>%
    mutate(covariate = fct_reorder(covariate, pre_weight)) %>% 
    pivot_longer(cols = c(pre_weight, post_weight), names_to = "balance", values_to = "mean_difference") %>% 
    mutate(allowability = ifelse(covariate %in% c(allowable_covs, "sexM", "sexF"), 1, 0)) %>% 
    mutate(covariate = fct_reorder(covariate, allowability)) %>% 
    slice(1:subset_covs) %>% 
    ggplot(aes(x = mean_difference, y = covariate, color = balance)) +
    geom_point(size = 6) +
    theme_minimal(base_size = 22) +
    scale_color_manual(values = c("pre_weight" = "red", "post_weight" = "blue")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(title = title, x = "Absolute Mean Difference", y = "", color = "Sample") +
    scale_x_continuous(limits = c(0, 1)) + # Set x-axis limits
    theme(legend.position = "bottom") + 
    theme(plot.title = element_text(hjust = 0.5)) #+ coord_flip()
}