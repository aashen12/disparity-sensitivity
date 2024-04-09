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