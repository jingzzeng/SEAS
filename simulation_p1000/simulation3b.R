# This file reproduces the results under model (M3b) in Table 1 with p = 1000.
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
# ------- Model 3b ------- #
p <- 1000    # Dimension of X
N <- 200     # Sample size
H <- 5       # Number of slices in SEAS-SIR and SEAS-Intra
s <- 6       # Sparsity level
d <- 2       # Structural dimension
ind <- 1:s   # Active set

load("beta/beta3b.rda")      # Load the beta matrix generated from "simulation3b_dat.R"

## main program
foo <- function(i){
  cat("Time", i, '\n')
  load(paste0("dat/sim3b/dat", i, ".rda"))        # Load the data set generated from "simulation3b_dat.R"
  load(paste0("dat/sim3b/foldid", i, ".rda"))     # Load the fold id generated from "simulation3b_dat.R"
  x_train <- dat$x
  y_train <- dat$y
  
  # ---------------- SEAS-SIR --------------- #
  seassir_fit <- cv.seas(x_train, y_train, H = H, foldid = foldid, type = 'sir')
  beta_seassir <- seassir_fit$beta
  if(!is.null(beta_seassir)){
    ind_hat <- which(apply(beta_seassir, 1, function(x) any(x!=0)))
    TPR_seassir <- sum(ind_hat %in% ind)/length(ind)
    FPR_seassir <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    d_seassir <- seassir_fit$rank
    dist_seassir <- subspace(beta, beta_seassir)
    seassir_result <- c(dist_seassir, TPR_seassir, FPR_seassir, d_seassir)
  }else{
    seassir_result <- rep(NA, 4)
  }
  names(seassir_result) <- c("dist_seassir", "TPR_seassir", "FPR_seassir", "d_seassir")
  
  # ---------------- SEAS-intra --------------- #
  seasintra_fit <- cv.seas(x_train, y_train, H = H, foldid = foldid, type = 'intra')
  beta_seasintra <- seasintra_fit$beta
  if(!is.null(beta_seasintra)){
    ind_hat <- which(apply(beta_seasintra, 1, function(x) any(x!=0)))
    TPR_seasintra <- sum(ind_hat %in% ind)/length(ind)
    FPR_seasintra <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    d_seasintra <- seasintra_fit$rank
    dist_seasintra <- subspace(beta, beta_seasintra)
    seasintra_result <- c(dist_seasintra, TPR_seasintra, FPR_seasintra, d_seasintra)
  }else{
    seasintra_result <- rep(NA, 4)
  }
  names(seasintra_result) <- c("dist_seasintra", "TPR_seasintra", "FPR_seasintra", "d_seasintra")
  
  # ---------------- SEAS-PFC --------------- #
  seaspfc_fit <- cv.seas(x_train, y_train, foldid = foldid, type = 'pfc')
  beta_seaspfc <- seaspfc_fit$beta
  if(!is.null(beta_seaspfc)){
    ind_hat <- which(apply(beta_seaspfc, 1, function(x) any(x!=0)))
    TPR_seaspfc <- sum(ind_hat %in% ind)/length(ind)
    FPR_seaspfc <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    d_seaspfc <- seaspfc_fit$rank
    dist_seaspfc <- subspace(beta, beta_seaspfc)
    seaspfc_result <- c(dist_seaspfc, TPR_seaspfc, FPR_seaspfc, d_seaspfc)
  }else{
    seaspfc_result <- rep(NA, 4)
  }
  names(seaspfc_result) <- c("dist_seaspfc", "TPR_seaspfc", "FPR_seaspfc", "d_seaspfc")
  
  c(seassir_result, seasintra_result, seaspfc_result)
}

times <- 100
output <- pbmclapply(seq_len(times), foo, mc.cores = 16)
# output <- lapply(seq_len(times), foo)
output <- do.call(rbind, output)
write.table(output, file = "output/M3b")
