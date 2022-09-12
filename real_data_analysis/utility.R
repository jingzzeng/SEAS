# Construct AR matrix
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# Cut small values in a matrix to zero. 
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}

# Evaluation based on distance correlation. (Szekely et al., 2007)
eval_dc <- function(Beta, x, y){
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      dcor(x %*% mat, y)
    }
  })
  return(result)
}

# Compute M and U matrices from observation data.
#########
# Input:
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in 1,...,H.
# type: Specifying the specific SEAS method. "sir" means SEAS-SIR, "intra" means SEAS-Intra and "pfc" means SEAS-PFC.
# FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3).
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
MU <- function(x, y, yclass=NULL, type='sir', FUN = NULL, categorical = FALSE, H = 5){
  if(is.null(yclass)){ # Construct the discretized response
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  cls <- sort(unique(yclass))
  nclass <- length(cls)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  mu <- colMeans(x)
  x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  M <- crossprod(x_c/sqrt(nobs)) # sample covariance of X
  
  if(type == 'sir'){
    U <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      U[, i] <- colMeans(x_c[yclass == cls[i],, drop=FALSE])
    }
  }else if(type == 'intra'){
    y_c <- y - mean(y)
    U <- matrix(0, nvars, nclass)
    lb <- quantile(y_c, 0.1)[[1]]
    ub <- quantile(y_c, 0.9)[[1]]
    y_c <- sapply(y_c, cut_func, lb = lb, ub = ub)
    for (i in 1:nclass){
      y_copy <- y_c
      y_copy[yclass!=cls[i]] <- 0
      U[, i] <- (1/nobs) * t(x_c) %*% (y_copy - mean(y_copy))
    }
  }else if(type == 'pfc'){
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) # the default function
    else Fmat <- t(sapply(y, FUN))
    Fmat_mean <- colMeans(Fmat)
    Fmat_c <- Fmat - matrix(Fmat_mean, NROW(Fmat), NCOL(Fmat), byrow = TRUE) # centered function f
    lb <- apply(Fmat_c, 2, quantile, 0.1)
    ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)){
      Fmat_c[,i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
    }
    U <- (1/nobs)*(t(x_c) %*% Fmat_c)
  }
  list(M = M, U = U, nclass = nclass, prior=prior)
}

# Cut extreme values in the samples. This function is used in MU function.
cut_func <- function(x, lb, ub){
  if(x < lb){
    return(lb)
  } else if(x > ub){
    return(ub)
  } else{
    return(x)
  }
}

# Estimate the rank of a matrix.
rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

# Subspace distance, defined in (19)
subspace <- function(A,B){
  if(is.vector(A)) A <- as.matrix(A)
  if(is.vector(B)) A <- as.matrix(B)
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

## ------------------------------------------------ ##
## The utility functions imported from R package 'msda'.
## These functions are used in 'msda' and 'cv.msda' functions.
formatoutput <- function(fit, maxit, pmax, p, H) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0){
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i){
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H)
      theta_i[ja,] <- a
      theta_i
    })
  }
  else{
    theta <- lapply(seq(nalam), function(x){matrix(0, p, H)})
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, " iterations; solutions for larger lambdas returned.\n", sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam){
  llam <- log(lam)
  if(length(llam) >= 3){lam[1] <- exp(2 * llam[2] - llam[3])}
  lam
}