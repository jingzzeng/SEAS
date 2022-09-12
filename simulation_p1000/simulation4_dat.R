# This file generate data from model (M4) with p=1000.
rm(list = ls())
library(pbmcapply)
library(MASS)
source("utility.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(2)

## Model setting
p <- 1000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 6       # Sparsity level
d <- 2       # Structural dimension
ind <- 1:s   # Active set

beta <- matrix(0, p, 2)                 # Mean of X
beta[ind,1] <- runif(s,2,2.5)
beta[c(1,3,5),2] <- runif(3,2,2.5)
beta[c(2,4,6),2] <- runif(3,-2.5,-2)
save(beta, file = "beta/beta4.rda")     # Save beta matrix.

Delta <- AR(0.5,p)                      # Covariance of the error term
Gamma <- Delta %*% beta                 # Represents $\Gamma \eta$ in Model (M4).

data_gen <- function(N){                # Data generating function
  y <- runif(N, -1, 1)
  f <- cbind(y, abs(y))
  error <- mvrnorm(N, rep(0,p), Delta)
  x <- f %*% t(Gamma) + error
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
  save(dat, file = paste0("dat/sim4/dat",i,'.rda'))     # Save data set.
  
  nfolds <- 5
  count <- rep(N/H, H)
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  save(foldid, file = paste0("dat/sim4/foldid",i,'.rda'))     # Save the cross-validation fold id.
}

times <- 100
pbmclapply(seq_len(times), foo, mc.cores = 16)
# lapply(1:times, foo)
