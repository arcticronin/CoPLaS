library(survival)
library(glmnet)
library(prioritylasso)

setwd("/home/ronin/Dev/notebooks/thesis_notebook/")

data <- read.csv("DaniDatasets/preprocessed.csv")
# remove first column with patient ID
data <- data[-c(1)] 


# Define the function
cox_lasso_model_with_offset <- function(df) {
  # Check for necessary columns
  required_cols <- c("time", "event")
  if (!all(required_cols %in% names(df))) {
    stop("Dataframe missing required columns: time, event")
  }
  
  # Prepare the survival object
  y <- with(df, Surv(time, event))
  
  # Prepare the matrix of predictors
  predictors <- setdiff(names(df), required_cols)
  x <- as.matrix(df[, predictors])
  
  # Extract the offset
  #offset_values <- df$offset
  
  # Fit Cox model with Lasso using cross-validation
  fit <- cv.glmnet(x, y, family="cox", alpha=1) 
  #fit <- cv.glmnet(x, y, family="cox", alpha=1, offset=offset_values) 
  
  # Best lambda and model fit
  best_lambda <- fit$lambda.min
  # best_fit <- glmnet(x, y, family="cox", alpha=1, lambda=best_lambda, offset=offset_values)
  best_fit <- glmnet(x, y, family="cox", alpha=1, lambda=best_lambda)
  
  # Extracting coefficients
  coef_vector <- as.vector(coef(best_fit, s = best_lambda))
  
  # handling intercept is not strictly recommended in COX, and not itnterpretable, 
  # so I'm leaving it out
  # coef_df <- data.frame(coefficient = coef_vector, row.names = c("(Intercept)", predictors))
  coef_df <- data.frame(coefficient = coef_vector, row.names = predictors)
  
  # Return results
  return(list(coef=coef_df, lambda=best_lambda, fit=best_fit))
}

model_results <- cox_lasso_model_with_offset(data)
 
#write.csv(model_results$coef, 
#          "results/glmnet_lasso_coeffs.csv", 
#          row.names = TRUE)


set.seed(123) # for reproducibility

folds <- cut(seq(1, nrow(data)), breaks=5, labels=FALSE) # 5-fold CV

# Store results
cv_results <- data.frame(fold = integer(0), c_index = numeric(0))

# start crossval and store C-Index
for(i in 1:5){
  # Split data into training and test sets
  test_indices <- which(folds == i)
  train_indices <- which(folds != i)
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit model on training data
  y_train <- with(train_data, Surv(time, event))
  x_train <- as.matrix(train_data[, setdiff(names(train_data), c("time", "event"))])
  fit <- cv.glmnet(x_train, y_train, family="cox", alpha=1)
  
  # Predict on test data
  x_test <- as.matrix(test_data[, setdiff(names(test_data), c("time", "event"))])
  predictions <- predict(fit, newx = x_test, s = "lambda.min")
  
  # Calculate C-index for test data
  c_index <- survConcordance(Surv(test_data$time, test_data$event) ~ predictions)$concordance
  
  # Save results
  cv_results <- rbind(cv_results, data.frame(fold = i, c_index = c_index))
}

# Calculate average C-index
mean_c_index <- mean(cv_results$c_index)



# Function to fit a Cox model with Ridge regularization on a range of columns
# used to score omics
cox_ridge_model <- function(df, start_col, end_col) {
  # Check for necessary columns as before
  required_cols <- c("time", "event")
  if (!all(required_cols %in% names(df))) {
    stop("Dataframe missing required columns: time, event")
  }
  
  # Validate column range
  if (start_col < 1 || end_col > length(df)) {
    stop("Invalid column range")
  }
  
  # Prepare the survival object
  y <- with(df, Surv(time, event))
  
  # Prepare the matrix of predictors
  predictors <- names(df)[start_col:end_col]
  x <- as.matrix(df[, predictors])
  
  # Fit Cox model with Ridge using cross-validation
  fit <- cv.glmnet(x, y, family="cox", alpha=0) # alpha=0 for Ridge
  
  # Best lambda and model fit
  best_lambda <- fit$lambda.min
  best_fit <- glmnet(x, y, family="cox", alpha=0, lambda=best_lambda)
  
  # Extracting coefficients
  coef_vector <- as.vector(coef(best_fit, s = best_lambda))
  coef_df <- data.frame(coefficient = coef_vector, row.names = predictors)
  
  # Return results
  return(list(coef=coef_df, lambda=best_lambda, fit=best_fit))
}

calc_omic_importance <- function(coef_df) {
  mean(sum(coef_df$coefficient^2))
}


indices_for_preprocessed <- read.csv("DaniDatasets/indices_for_preprocessed.csv")


automate_ridge_fitting <- function(data, indices_for_preprocessed) {
  indices <- as.numeric(indices_for_preprocessed$index)
  
  scores <- list()
  
  for (i in 1:(length(indices) - 1)) {
    start_col <- indices[i] + 1
    end_col <- indices[i + 1]
    
    # Fit the Ridge model with columns
    ridge_fit <- cox_ridge_model(data, start_col, end_col)
    
    score <- calc_omic_importance(ridge_fit$coef)
    scores[[i]] <- score
  }
  
  names(scores) <- c("clinical", 
                     "rna", 
                     "methylation", 
                     "mirna", 
                     "cna_log2", 
                     "rppa", 
                     "alterations", 
                     "microbiome")
  return(scores)
}

# Usage
scores <- automate_ridge_fitting(data, indices_for_preprocessed)
