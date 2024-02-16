library(survival)
library(glmnet)

setwd("/home/ronin/Dev/notebooks/thesis_notebook/")
# Assuming you have already set your working directory and loaded the data
step2 <- read.csv("step2.csv")

# Define the function
cox_lasso_model_with_offset <- function(df) {
  # Check for necessary columns
  required_cols <- c("time", "event", "offset")
  if (!all(required_cols %in% names(df))) {
    stop("Dataframe missing required columns: time, event, or offset")
  }
  
  # Prepare the survival object
  y <- with(df, Surv(time, event))
  
  # Prepare the matrix of predictors
  predictors <- setdiff(names(df), required_cols)
  x <- as.matrix(df[, predictors])
  
  # Extract the offset
  offset_values <- df$offset
  
  # Fit Cox model with Lasso using cross-validation
  fit <- cv.glmnet(x, y, family="cox", alpha=1, offset=offset_values) 
  
  # Best lambda and model fit
  best_lambda <- fit$lambda.min
  best_fit <- glmnet(x, y, family="cox", alpha=1, lambda=best_lambda, offset=offset_values)
  
  # Extracting coefficients
  coef_vector <- as.vector(coef(best_fit, s = best_lambda))
  
  # handling intercept is not strictly recommended in COX, and not itnterpretable, 
  # so I'm leaving it out
  # coef_df <- data.frame(coefficient = coef_vector, row.names = c("(Intercept)", predictors))
  coef_df <- data.frame(coefficient = coef_vector, row.names = predictors)
  
  # Return results
  return(list(coef=coef_df, lambda=best_lambda, fit=best_fit))
}

# Example usage:
model_results <- cox_lasso_model_with_offset(step2)
head(model_results$coef)
