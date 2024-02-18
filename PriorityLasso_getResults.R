# survival
library(survival)
library(glmnet)
library(prioritylasso)

# plotting
library(ggplot2)
library(reshape2)

# confidence interval
library(dplyr)
library(tidyr)

# todo, aggiungi il preprocessing qui
 

set.seed(42)

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
  
  #fit <- TODO aggiungi il fit qui
  
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

