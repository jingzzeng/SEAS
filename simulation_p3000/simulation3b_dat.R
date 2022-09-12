# This file generate data from model (M3b) with p=3000.
rm(list = ls())
library(pbmcapply)
library(MASS)
source("utility.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(3)

## Model setting
p <- 3000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 6       # Sparsity level
d <- 2       # Structural dimension
ind <- 1:s   # Active set

Mu <- list()      # Mean vectors for Gaussian mixture model.
Mu1 <- rep(0, p); Mu1[ind] <- rep(-1, s)
Mu2 <- rep(0, p); Mu2[ind] <- rep(0, s)
Mu3 <- rep(0, p); Mu3[ind] <- rep(1, s)
Mu[[1]] <- Mu1; Mu[[2]] <- Mu2; Mu[[3]] <- Mu3
Sigma <- list()   # Covariance matrices for Gaussian mixture model.
Sigma[[1]] <- AR(0.1, p)
Sigma[[2]] <- AR(0.5, p)
Sigma[[3]] <- AR(0.9, p)
beta <- matrix(0, p, 2)     # Central subspace basis
beta[ind,1] <- runif(s,0.3, 0.6)
beta[ind,2] <- c(runif(3, 0.3, 0.6), runif(3,-0.6,-0.3))
save(beta, file = "beta/beta3b.rda")    # Save beta matrix.

data_gen <- function(N){    # Data generating function
  I <- sample(1:3, N, replace = TRUE, prob = c(0.4, 0.2, 0.4))
  x <- matrix(0, N, p)
  for(i in 1:3){
    x[I == i,] <- mvrnorm(sum(I == i), Mu[[i]], Sigma[[i]])
  }
  y <- x %*% beta[,1] * exp(x %*% beta[,2] + 0.5 *  rnorm(N))
  list(x = x, y = y)
}

## ------------ Main function ------------ ##
foo <- function(i){
  cat("Time", i, '\n')
  train <- data_gen(N)          # Generate data set
  x_train <- train$x
  y_train <- train$y  
  ord <- order(y_train)         # Sort data set
  y_train <- y_train[ord]
  x_train <- x_train[ord,]
  
  dat <- list(x = x_train, y = y_train)
  save(dat, file = paste0("dat/sim3b/dat",i,'.rda'))    # Save data set.
  
  nfolds <- 5
  count <- rep(N/H, H)
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  save(foldid, file = paste0("dat/sim3b/foldid",i,'.rda'))     # Save the cross-validation fold id.
}

times <- 100
pbmclapply(seq_len(times), foo, mc.cores = 16)
# lapply(1:times, foo)
