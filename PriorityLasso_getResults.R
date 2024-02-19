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


setwd("Dev/notebooks/thesis_notebook/DaniDatasets/")
data = read.csv("preprocessed.csv")
indices = read.csv("indices_for_preprocessed.csv")

#preproc, remove dependent var and iupdate indices
X <- as.matrix(data[, !(names(data) %in% c("time", "event", "X"))])  # Predictor matrix
y <- Surv(data$time, data$event)  # Response variable for survival analysis

indices$index[2:nrow(indices)] <- indices$index[2:nrow(indices)] - 2

blocks = list(bp1=1:2, 
              bp2=3:1999, 
              bp3=2000:2487, 
              bp4=2488:3188, 
              bp5=3189:5464, 
              bp6=5465:5587, 
              bp7=5588:7036, 
              bp8=7037:8376)

# Creare una lista per memorizzare le liste di indici
list_of_indices <- list()

set.seed(42)

folds <- cut(seq(1, nrow(data)), breaks=30, labels=FALSE) # 5-fold CV

# Store results
cv_results <- data.frame(fold = integer(0), c_index = numeric(0))

# start crossval and store C-Index
for(i in 1:30){
  # Split data into training and test sets
  test_indices <- which(folds == i)
  train_indices <- which(folds != i)
  train_data <- X[train_indices, ]
  test_data <- X[test_indices, ]
  
  # Fit model on training data
  y_train = y[train_indices,]
  x_train = X[train_indices,]
  
  fit <- prioritylasso(Y = y_train, 
                       X = x_train, 
                       blocks = blocks, 
                       family = "cox", 
                       block1.penalization = TRUE, 
                       type.measure = "deviance",
                       lambda.type = "lambda.min",
                       nfolds = 10)
  
  # Predict on test data
  x_test <- X[test_indices,]
  y_test = y[test_indices,]
  
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = "all"
  )
  
  # print(predictions)
  # print(y[test_indices,])
  
  # Calculate C-index for test data
  # c_index <- concordance(y_test ~ predictions)$concordance
  #c_index <- concordance(y_test, predictions)
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results <- rbind(cv_results, data.frame(fold = i, c_index = c_index))
}

# Calculate average C-index
mean_c_index <- mean(cv_results$c_index)
print(mean_c_index)

