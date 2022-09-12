### Lymphoma data description ###
# The lymphoma dataset consists of 42 samples of diffuse large B-cell lymphoma (DLBCL), 9 samples of follicular lymphoma (FL), and 11 samples of chronic lymphocytic leukemia (CLL). DBLCL, FL, and CLL classes are coded in 0, 1, and 2, respectively, in y vector. Matrix x is gene expression data and arrays were normalized, imputed, log transformed, and standardized to zero mean and unit variance across genes.
####################
# This file reproduces the estimation results in Table 2.
rm(list = ls())
library(pbmcapply)
library(LassoSIR)
library(energy)
library(glmnet)
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

######## Estimation consistency ##########
foo <- function(x, y){
  seassir_fit <- cv.seas(x, y, categorical=TRUE, foldid = foldid, lam1_fac=seq(1.2,0.5, length.out = 10), type = 'sir')
  beta_seassir <- seassir_fit$beta
  if(is.vector(beta_seassir)) beta_seassir <- matrix(beta_seassir)
  if(is.null(beta_seassir) || (seassir_fit$rank==0)){
    d_seassir <- NA
    ind_hat_seassir <- NA
    s_seassir <- NA
  }else{
    d_seassir <- seassir_fit$rank
    ind_hat_seassir <- which(apply(beta_seassir, 1, function(x) any(x!=0)))
    s_seassir <- length(ind_hat_seassir)
  }
  
  lassosir_fit <- LassoSIR_revised(x, y, categorical = TRUE, foldid = foldid, choosing.d = 'automatic')
  beta_lassosir <- lassosir_fit$beta
  if(is.vector(beta_lassosir)) beta_lassosir <- matrix(beta_lassosir)
  if(is.null(beta_lassosir) || (lassosir_fit$no.dim==0)){
    d_lassosir <- NA
    ind_hat_lassosir <- NA
    s_lassosir <- NA
  }else{
    d_lassosir <- lassosir_fit$no.dim
    ind_hat_lassosir <- which(apply(beta_lassosir, 1, function(x) any(x!=0)))
    s_lassosir <- length(ind_hat_lassosir)
  }
  
  output <- list(rank = c(d_seassir, d_lassosir), s = c(s_seassir, s_lassosir), ind_hat = list(ind_hat_seassir, ind_hat_lassosir), beta = list(beta_seassir, beta_lassosir))
  output
}

# Bootstrap function
boot_foo <- function(i){
  cat('Time', i, '\n')
  # Boostrap samples.
  index <- sort(sample(1:length(y), length(y), replace = TRUE))
  boot_x <- x[index,]
  boot_y <- y[index]
  # Re-sample if any counts less than nfold = 5
  while(any(table(boot_y) < 5)){
    index <- sort(sample(1:length(y), length(y), replace = TRUE))
    boot_x <- x[index,]
    boot_y <- y[index]
  }
  
  # Fold id
  nfolds <- 5
  count <- as.numeric(table(boot_y))
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  }
  
  # Estimate the central subspace using the bootstrap samples.
  boot_output <- foo(boot_x, boot_y)
  rank <- boot_output$rank
  s <- boot_output$s
  ind_hat <- boot_output$ind_hat
  beta <- boot_output$beta
  # Compute the simple matching coefficient.
  smc <- sapply(1:length(ind_hat), function(k){
    if (is.na(ind_hat[[k]]) || is.na(output$ind_hat[[k]])) 0
    else 1-(length(union(output$ind_hat[[k]], ind_hat[[k]])) - length(intersect(output$ind_hat[[k]], ind_hat[[k]])))/NCOL(x)
  })
  # Compute the subspace distance between the subspaces estimated from both the original data and the bootstrap samples.
  dist <- sapply(1:length(beta), function(k){
    if(is.null(output$beta[[k]]) || is.null(beta[[k]])) NA
    else subspace(output$beta[[k]], beta[[k]])
  })
  out <- paste0(dist, " ", smc, " ", s, " ", rank, collapse = " ")
  as.numeric(strsplit(out, " ")[[1]])
}

times <- 100
nfolds <- 5
count <- as.numeric(table(y))
foldid <- c()
for(cnt in count){
  foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
}
# Estimate the central subspace from the original data set.
output <- foo(x,y)
# Estimate the central subspace from the bootstrap samples.
boot_output <- pbmclapply(1:times, boot_foo, mc.cores = 16)
# boot_output <- lapply(seq_len(times), boot_foo)
boot_output <- do.call(rbind, boot_output)
write.table(boot_output, file = "output/lymphoma_est")
