# This file reproduces the results for Lasso-SIR under model (M2) in Table 1 with p = 1000.
rm(list = ls())
library(pbmcapply)
library(MASS)
library(glmnet)
library(LassoSIR)
library(ggplot2)
library(energy)
library(msda)
source("utility.R")
source("seas.R")
source("LassoSIR_revised.R")
# #############################################
# ------- Model 2 ------- #
p <- 1000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 10      # Sparsity level
d <- 1       # Structural dimension
ind <- 1:s   # Active set

load("beta/beta2.rda")         # Load the beta matrix generated from "simulation2_dat.R"

## main program
foo <- function(i){
  cat("Time", i, '\n')
  load(paste0("dat/sim2/dat", i, ".rda"))        # Load the data set generated from "simulation2_dat.R"
  load(paste0("dat/sim2/foldid", i, ".rda"))     # Load the fold id generated from "simulation2_dat.R"
  x_train <- dat$x
  y_train <- dat$y
  
  # ---------------- LassoSIR --------------- #
  lassosir_fit <- LassoSIR_revised(x_train, y_train, H = H, foldid = foldid, choosing.d = 'automatic')
  beta_lassosir <- lassosir_fit$beta
  if(!is.null(beta_lassosir)){
    ind_hat <- which(apply(beta_lassosir, 1, function(x) any(x!=0)))
    TPR_lassosir <- sum(ind_hat %in% ind)/length(ind)
    FPR_lassosir <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    d_lassosir <- lassosir_fit$no.dim
    dist_lassosir <- subspace(beta, beta_lassosir)
    lassosir_result <- c(dist_lassosir, TPR_lassosir, FPR_lassosir, d_lassosir)
  }else{
    lassosir_result <- rep(NA, 4)
  }
  names(lassosir_result) <- c("dist_lassosir", "TPR_lassosir", "FPR_lassosir", "d_lassosir")
  
  lassosir_result
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
times <- 100
output <- pbmclapply(seq_len(times), foo, mc.cores = 16)
# output <- lapply(seq_len(times), foo)
output <- do.call(rbind, output)
write.table(output, file = "output/M2_LassoSIR")
