# This file generate data from model (M2) with p=3000.
rm(list = ls())
library(pbmcapply)
library(MASS)
source("utility.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(2)

## Model setting
p <- 3000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 10                           # Sparsity level
d <- 1                            # Structural dimension
ind <- 1:s                        # Active set

Mu <- rep(0,p)                    # Mean of X
Sigma <- AR(0.5, p)               # Covariance of X
beta <- matrix(0, p, 1)           # Central subspace basis
beta[ind,1] <- runif(s, 0.1, 0.4) 
save(beta, file = "beta/beta2.rda")  # Save beta matrix.

data_gen <- function(N){          # Data generating function
  x <- mvrnorm(N, Mu, Sigma)
  nobs <- dim(x)[1]
  y <- sinh(x %*% beta) + rnorm(nobs)
  list(x = x, y = y)
}

## ------------ Main function ------------ ##
foo <- function(i){
  cat("Time", i, '\n')
  train <- data_gen(N)         # Generate data set
  x_train <- train$x
  y_train <- train$y  
  ord <- order(y_train)        # Sort data set
  y_train <- y_train[ord]
  x_train <- x_train[ord,]
  
  dat <- list(x = x_train, y = y_train)
  save(dat, file = paste0("dat/sim2/dat",i,'.rda'))    # Save data set.
  
  nfolds <- 5
  count <- rep(N/H, H)
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  save(foldid, file = paste0("dat/sim2/foldid",i,'.rda'))    # Save the cross-validation fold id.
}

times <- 100
pbmclapply(seq_len(times), foo, mc.cores = 16)
# lapply(1:times, foo)
