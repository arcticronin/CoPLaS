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


#setwd("Dev/notebooks/thesis_notebook/DaniDatasets/")
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
cv_results_1 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_12 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_13 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_123 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_1367 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_13678 <- data.frame(fold = integer(0), c_index = numeric(0))
cv_results_3678 <- data.frame(fold = integer(0), c_index = numeric(0))

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
  
  # 1
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_1 <- rbind(cv_results_1, data.frame(fold = i, c_index = c_index))
  
  # 1 e 2
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1, 2))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_12 <- rbind(cv_results_12, data.frame(fold = i, c_index = c_index))
  
  # 1 e 3
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1, 3))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_13 <- rbind(cv_results_13, data.frame(fold = i, c_index = c_index))
  
  # 1 2 e 3
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1, 2, 3))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_123 <- rbind(cv_results_123, data.frame(fold = i, c_index = c_index))
  
  # 1 3 6 7
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1, 3, 6, 7))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_1367 <- rbind(cv_results_1367, data.frame(fold = i, c_index = c_index))
  
  # 1 3 6 7 8
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(1, 3, 6, 7, 8))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_13678 <- rbind(cv_results_13678, data.frame(fold = i, c_index = c_index))
  
  # 3 6 7 8
  predictions <- predict(
    object = fit, 
    newdata = x_test, 
    type = "response",
    use.blocks = c(3, 6, 7, 8))
  c_index <- survConcordance(y_test ~ predictions)$concordance
  # Save results
  cv_results_3678 <- rbind(cv_results_3678, data.frame(fold = i, c_index = c_index))
}

# Calculate average C-index
cv_res <- cv_results_1[complete.cases(cv_results_1), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.600
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_12[complete.cases(cv_results_12), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.62
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_13[complete.cases(cv_results_13), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.59
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_123[complete.cases(cv_results_123), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.612
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_1367[complete.cases(cv_results_1367), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.614
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_13678[complete.cases(cv_results_13678), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.63
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]

cv_res <- cv_results_3678[complete.cases(cv_results_3678), ]
cv_res <- cv_res[rowSums(cv_res == 0) == 0, ]
mean_c_index <- mean(cv_res$c_index)
print(mean_c_index) #0.5959
ci <- t.test(cv_res$c_index, conf.level = 0.95)
ci$conf.int[1]
ci$conf.int[2]
