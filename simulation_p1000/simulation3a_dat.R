# This file generate data from model (M3a) with p=1000.
rm(list = ls())
library(pbmcapply)
library(MASS)
source("utility.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

## Model setting
p <- 1000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 6       # Sparsity level
d <- 2       # Structural dimension
ind <- 1:s   # Active set

Mu <- rep(0,p)                    # Mean of X
Sigma <- AR(0.5, p)               # Covariance of X
beta <- matrix(0, p, 1)           # Central subspace basis
beta[ind,1] <- runif(s,0.3, 0.6)
beta[ind,2] <- c(runif(3, 0.3, 0.6), runif(3,-0.6,-0.3))
save(beta, file = "beta/beta3a.rda")   # Save beta matrix.

data_gen <- function(N){          # Data generating function
  x <- mvrnorm(N, Mu, Sigma)
  nobs <- dim(x)[1]
  y <- x %*% beta[,1] * exp(x %*% beta[,2] + 0.5 * rnorm(nobs))
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
  save(dat, file = paste0("dat/sim3a/dat",i,'.rda'))    # Save data set.
  
  nfolds <- 5
  count <- rep(N/H, H)
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  save(foldid, file = paste0("dat/sim3a/foldid",i,'.rda'))     # Save the cross-validation fold id.
}

times <- 100
pbmclapply(seq_len(times), foo, mc.cores = 16)
# lapply(1:times, foo)
