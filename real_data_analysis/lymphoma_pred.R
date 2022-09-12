### Lymphoma data description ###
# The lymphoma dataset consists of 42 samples of diffuse large B-cell lymphoma (DLBCL), 9 samples of follicular lymphoma (FL), and 11 samples of chronic lymphocytic leukemia (CLL). DBLCL, FL, and CLL classes are coded in 0, 1, and 2, respectively, in y vector. Matrix x is gene expression data and arrays were normalized, imputed, log transformed, and standardized to zero mean and unit variance across genes.
####################
# This file reproduces the classification results in Table 2.
rm(list = ls())
library(pbmcapply)
library(MASS)
library(LassoSIR)
library(nnet)
library(glmnet)
library(e1071)
library(randomForest)
library(energy)
library(R.matlab)
library(msda)
source("utility.R")
source("seas.R")
source("LassoSIR_revised.R")
data("lymphoma", package = "spls")
#############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

x <- lymphoma$x
y <- lymphoma$y
y <- y + 1

# Construct dummy variable for y.
y_dummy <- sapply(unique(y), function(i){
  as.numeric(y == i)
})

# Variable screen: keep 200 variables.
dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
  dcor(y_dummy, x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:200]
x <- x[,ord, drop=FALSE]
ratio <- 0.8

# Prediction function
err_func <- function(x_train, y_train, x_test, y_test){
  x_train <- data.frame(x_train)
  x_test <- data.frame(x_test)
  colnames(x_train) <- paste0('X', seq_len(ncol(x_train)))
  colnames(x_test) <- paste0('X', seq_len(ncol(x_test)))
  data <- cbind(y = factor(y_train), x_train)
  
  # Multinomial logistic regression
  model <- nnet::multinom(y~., data = data, trace = FALSE)
  prediction <- predict(model, newdata=x_test, type='class')
  log_err <- mean(prediction != y_test)
  
  # Support vector machine.
  model <- e1071::svm(y ~., data = data, scale = FALSE, kernel = 'linear', cost = 1, probability = TRUE)
  prediction <- predict(model, newdata = x_test, probability = TRUE)
  svm_err <- mean(prediction != y_test)
  
  # Linear discriminant analysis.
  model <- tryCatch(MASS::lda(y~., data = data),
                    error = function(e){
                      return(NULL)
                    })
  if(is.null(model)){
    prediction <- rep(unique(y_test)[which.max(table(y_test))], length(y_test))
  }
  else prediction <- predict(model, newdata = x_test)$class
  lda_err <- mean(prediction != y_test)
  
  # Random forest
  model <- randomForest::randomForest(y~., data = data)
  prediction <- predict(model, newdata = x_test)
  rf_err <- mean(prediction != y_test)
  
  out <- c(log_err = log_err, svm_err = svm_err, lda_err = lda_err, rf_err = rf_err)
}

# Main function
foo <- function(i){
  cat(c('Time', i, '\n'))
  class <- unique(y)
  index <- c()
  ## Stratified folding
  for(k in class){
    index <- c(index, sample(which(y==k), ratio*length(y[y==k]), replace = FALSE))
  }
  
  x_train <- x[index,,drop=FALSE]
  y_train <- y[index]
  x_test <- x[-index,,drop=FALSE]
  y_test <- y[-index]
  
  # Fold id
  nfolds <- 5
  count <- as.numeric(table(y_train))
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  
  seassir_fit <- cv.seas(x_train, y_train, categorical=TRUE, foldid = foldid, lam1_fac=seq(1.2,0.5, length.out = 10), type = 'sir')
  beta_seassir <- seassir_fit$beta
  if(is.null(beta_seassir) || (seassir_fit$rank==0)){
    err_seassir <- rep(NA, 4)
  }else{
    x_train_new <- x_train %*% beta_seassir                              # Reduce the training set.
    x_test_new <- x_test %*% beta_seassir                                # Reduce the test set.
    err_seassir <- err_func(x_train_new, y_train, x_test_new, y_test)    # The classification error.
  }
  names(err_seassir) <- c("log_seassir", "svm_seassir", "lda_seassir", "rf_seassir")
  
  lassosir_fit <- LassoSIR_revised(x_train, y_train, categorical = TRUE, foldid = foldid, choosing.d = 'automatic')
  beta_lassosir <- lassosir_fit$beta
  if(is.null(beta_lassosir) || (lassosir_fit$no.dim==0)){
    err_lassosir <- rep(NA, 4)
  }else{
    x_train_new <- x_train %*% beta_lassosir
    x_test_new <- x_test %*% beta_lassosir
    err_lassosir <- err_func(x_train_new, y_train, x_test_new, y_test)
  }
  names(err_lassosir) <- c("log_lassosir", "svm_lassosir", "lda_lassosir", "rf_lassosir")
  
  c(err_seassir, err_lassosir)
}

times <- 100
output <- pbmclapply(seq_len(times), foo, mc.cores = 16)
# output <- lapply(seq_len(times), foo)
output <- do.call(rbind, output)
write.table(output, file = "output/lymphoma_pred")
