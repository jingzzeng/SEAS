## ---------------------------------------- ##
# This is the main file containing the functions to implement SEAS algorithm (Algorithm 1):
# #####
# seas: Estimate the central subspace using SEAS algorithm.
# cv.seas: The cross-validation function for SEAS algorithm.
# admm: The ADMM algorithm for solving Eq. (13) in Section 3.1.
# updataB: Update B matrix in Step (3a) of ADMM algorithm.
# updateC: Update C matrix in Step (3b) of ADMM algorithm.
# msda, cv.msda: The revised 'msda' and 'cv.msda' functions from package 'msda'. These two functions are used for automatically generating tuning parameter sequences in 'seas' and 'cv.seas' functions.
## ---------------------------------------- ##

seas <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', M = NULL, U = NULL, nobs = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), FUN = NULL, eps = 1e-3, maxit = 1e+3, ...){
  # Inputs:
  # =======
  # x: n x p observation matrix for predictor.
  # y: n-dimensional observation vector for response.
  # yclass: Discretized response taking values in 1,...,H.
  # d: True structural dimension. The default is NULL.
  # categorical: A logical value indicating whether y is categorical.
  # H: The number of slices. The default value is 5.
  # type: Specifying the specific SEAS method. "sir" means SEAS-SIR, "intra" means SEAS-Intra and "pfc" means SEAS-PFC.
  # M: The M matrix in optimization problem. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC. If both arguments 'M' and 'U' are provided, the arguments 'x' and 'y' are ignored.
  # U: The U matrix in optimization problem. If both arguments 'M' and 'U' are provided, the arguments 'x' and 'y' are ignored.
  # nobs: The number of observations.
  # lam1, lam2, gamma: The user-provided sequences of tuning parameters lambda_1, lambda_2 and gamma. If NULL, the sequences are automatically generated. 
  # lam1_fac, lam2_fac: The factors used in automatically generating the tuning parameter sequences.
  # FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3).
  # eps: The tolerance of convergence in ADMM algorithm. The value is passed to 'admm' function.
  # maxit: The maximal iterations in ADMM algorithm. The value is passed to 'admm' function.
  # 
  # Outputs:
  # ========
  # beta: A list containing the estimated basis matrices of central subspace.
  # B: A list containing the estimated B matrices.
  # rank: A vector containing the estimated ranks.
  # s: A vector containing the estimated sparsity levels.
  # step: A vector containing the number of iterations to converge for each tuning parameter.
  # lam1, lam2, gamma: The tuning parameter sequences.
  # code: Error code. If code == 0, exit normally.
  
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
        nclass <- as.integer(length(unique(yclass)))
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    if(is.null(gamma)){
      gamma <- c(10,30,50)
    }
    if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
      fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=10, lambda.factor=0.5, nfolds = 5, FUN = FUN, maxit=1e3)
      M <- fit_1$M  # The M matrix based on the full data.
      U <- fit_1$U  # The U matrix based on the full data.
      id_max_msda <- fit_1$id
      lam1_max_msda <- fit_1$lam_max  # The optimal lambda from msda
      beta_msda <- as.matrix(fit_1$beta)  # The optimal matrix from msda
      if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
      if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
      if (all(lam2 == 0)){
        lam2 <- 0
        warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
      }
    }else{
      MU_out <- MU(x, y, yclass, type, FUN)
      M <- MU_out$M
      U <- MU_out$U
    }
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(lam1) || is.null(lam2) || is.null(gamma)) stop("Sequences lam1, lam2 or gamma is missing.")
    if(is.null(nobs)) stop("Missing nobs.")
    nvars <- NCOL(M)
  }
  
  ## Error code
  code <- 0
  
  if(is.vector(lam1) && (length(lam1) == 1) && (lam1 == 0) && is.vector(lam2) && (length(lam2) == 1) && (lam2 == 0)){ # For degenerate case where lambda1 = lambda2 = 0, return B = M^{-1} U directly.
    B <- solve(M) %*% U
    if(is.null(d)) beta <- svd(B)$u
    else if(d == 0) beta <- matrix(0, nrow(Bnew), ncol(Bnew))
    else beta <- svd(B)$u[,1:d,drop=FALSE]
    vec <- as.vector(beta)
    vec[abs(vec) < 1e-3] <- 0
    beta <- matrix(vec, nrow(beta), ncol(beta))
    rank <- NCOL(beta)
    output <- list(beta = beta, B = B, rank = rank, lam1 = lam1, lam2 = lam2, code = code)
  }
  else{
    # Fit with admm function
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, ...)
    B_l <- fit$B
    beta_l <- fit$beta
    if (all(sapply(beta_l, is.null))){
      code <- 1
      warning("No converged results returned.")
      return(list(beta = beta_l, code = code))
    }
    rank_l <- fit$rank
    s_l <- fit$s
    step_l <- fit$step
    time_l <- fit$time
    if(length(B_l) == 1){
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}

cv.seas <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', lambda.factor=0.5, nlambda=10, nfolds = 5, foldid = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), plot = FALSE, FUN = NULL, eps = 1e-3, maxit = 1e+3, trace.it = FALSE, ...){
  # The inputs and outputs are similar to the ones in 'seas' functions. Only the different ones are listed below.
  # Inputs:
  # =======
  # nfolds: The number of folds in the cross-validation.
  # plot: If TRUE, (1) plot the evaluation for each tuning parameter in 'msda' function; (2) in each cross-validation data fold, plot the evaluation for each tuning parameter in 'seas' function.
  # trace.it: If TRUE, print the process of cross-validation.
  # 
  # Outputs:
  # ========
  # Refer to the outputs in 'seas' function.
  
  start_time <- Sys.time()
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  
  M <- U <- M_fold <- U_fold <- NULL
  nobs <- dim(x)[1]
  
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    if (nfolds < 3) stop("nfolds must be larger than 3")
    if (nfolds > nobs) stop("nfolds is larger than the sample size")
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(foldid))
  }
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
    fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=nlambda, lambda.factor=lambda.factor, foldid = foldid, FUN = FUN, maxit=1e3, plot = plot)
    M <- fit_1$M  # The M matrix based on the full data.
    U <- fit_1$U  # The U matrix based on the full data.
    M_fold <- fit_1$M_fold  # The M matrix based on each four out of five folds.
    U_fold <- fit_1$U_fold  # The U matrix based on each four out of five folds.
    id_max_msda <- fit_1$id
    lam1_max_msda <- fit_1$lam_max
    beta_msda <- as.matrix(fit_1$beta)
    if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)){
      lam2 <- 0
      warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # The number of errors
  nerr <- 0
  code <- 0
  
  end_time <- Sys.time()
  time1 <- difftime(end_time, start_time, units = "secs")
  ## Record time1: estimate M and U matrices
  
  out_all <- lapply(1:nfolds, function(k){ # Cross-validation
    if(trace.it) cat(sprintf("Fold: %d/%d\n", k, nfolds))
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]

    if(is.null(M_fold) || is.null(U_fold)){
      x_train <- x[foldid!=k,,drop=FALSE]
      y_train <- y[foldid!=k]
      yclass_train <- yclass[foldid!=k]
      # Fit with seas function
      fit_fold <- seas(x_train, y_train, yclass = yclass_train, type = type, FUN = FUN, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
    }
    else{
      fit_fold <- seas(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(foldid!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
    }

    err <- 0
    if(fit_fold$code != 0){
      err <- 1
      warning(paste0("Fold ", k, ": No converged results returned from the seas function."))
      eval_fold <- rep(NA_real_, length(fit_fold$beta))
      out <- list(eval_fold, err)
      return(out)
    }

    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time

    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)

    if(plot){ # If true, plot the evaluation for each tuning parameter
      dat <- data.frame(x = 1:length(eval_fold), y = eval_fold, rank = as.factor(rank_l))
      g <- ggplot(dat, aes(x = x, y = y, col = rank))+
        geom_point(size = 1)+
        labs(title=paste0("Fold ", k), x="", y="Distance correlation")+
        theme_bw()
      print(g)
    }
    out <- list(eval_fold, err)
    out
  })
  
  # Combine the evaluations from each fold.
  eval_all <- do.call(rbind, lapply(out_all, "[[", 1))
  errs <- do.call(c, lapply(out_all, "[[", 2))
  nerr <- sum(errs)
  
  if((nerr != 0) && (nerr != nfolds)){
    code <- 3
    warning(paste0("No converged results returned in", nerr, "folds."))
  }
  else if(nerr == nfolds){
    code <- 4
    warning("No converged results returned in any fold.")
    return(list(beta = NULL, code = code))
  }
  
  if(is.vector(eval_all)){
    eval_all <- as.matrix(eval_all)
  }
  
  # Compute the cross-validation mean and standard error.
  cvm <- colMeans(eval_all, na.rm=TRUE)
  cvsd <- sqrt(colMeans(scale(eval_all, cvm, FALSE)^2, na.rm = TRUE)/(nfolds-1))
  
  # Select the best tuning parameter.
  id_max <- which.max(cvm)
  id_lam1 <- ceiling(id_max/(n2*n3))
  id_gamma <- ceiling((id_max-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_max-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_max <- lam1[id_lam1]
  gamma_max <- gamma[id_gamma]
  lam2_max <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  start_time <- Sys.time()
  # Refit with the selected tuning parameters.
  if(is.null(M) || is.null(U)){
    fit <- seas(x, y, yclass = yclass, type = type, FUN = FUN, lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }
  else{
    fit <- seas(M = M, U = U, nobs = NROW(x), lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }
  
  if(fit$code != 0){
    code <- 5
    warning("The estimated beta is null.")
    return(list(beta = NULL, code = code))
  }
  
  B <- fit$B
  beta <- fit$beta
  rank <- fit$rank
  
  end_time <- Sys.time()
  time2 <- difftime(end_time, start_time, units = "secs") # We do not include the time for tuning parameter selection.
  # Record time: one run with the selected tuning parameter
  
  time <- time1 + time2
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval_all, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_max = lam1_max, lam2_max = lam2_max, gamma_max = gamma_max, code = code, time = time)
  output
}

# admm algorithm function
admm <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, ...){
  # Inputs:
  # =======
  # M: The M matrix in optimization problem. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
  # U: The U matrix in optimization problem.
  # nobs: The number of observations.
  # nvars: The number of predictors.
  # lam1, lam2, gamma: The user-specified sequences of tuning parameter lambda_1, lambda_2 and gamma.
  # eps: The tolerance of convergence in ADMM algorithm. The value is passed to 'admm' function.
  # maxit: The maximal iterations in ADMM algorithm. The value is passed to 'admm' function.
  # d: The true structural dimension. The default is NULL.
  # 
  # Outputs:
  # ========
  # beta: A list containing the estimated basis matrices of central subspace.
  # B: A list containing estimated B matrices.
  # rank: A vector containing the estimated ranks.
  # s: A vector containing the estimated sparsity levels.
  # step: A vector containing the number of iterations to converge for each tuning parameter.
  # nlam: The number of converged matrices.
  
  # since the user is required to provide lam1, then set flmin=1
  opts <- list(...)
  if(is.null(opts$nlam)) opts$nlam <- as.integer(1)
  if(is.null(opts$H)) opts$H <- as.integer(dim(U)[2])
  if(is.null(opts$nvars)) opts$nvars <- as.integer(nvars)
  if(is.null(opts$pf)) opts$pf <- as.double(rep(1, nvars))
  if(is.null(opts$dfmax)) opts$dfmax <- as.integer(nobs)
  if(is.null(opts$pmax)) opts$pmax <- as.integer(min(nobs*2+20, nvars))
  if(is.null(opts$flmin)) opts$flmin <- as.double(1)
  if(is.null(opts$eps_inner)) opts$eps_inner <- as.double(1e-04)
  if(is.null(opts$maxit_inner)) opts$maxit_inner <- as.integer(1e+6)
  if(is.null(opts$sml)) opts$sml <- as.double(1e-6)
  if(is.null(opts$verbose)) opts$verbose <- as.integer(FALSE)
  if(is.null(opts$nalam)) opts$nalam <- integer(1)
  if(is.null(opts$theta)) opts$theta <- double(opts$pmax * opts$H * opts$nlam)
  if(is.null(opts$itheta)) opts$itheta <- integer(opts$pmax)
  if(is.null(opts$ntheta)) opts$ntheta <- integer(opts$nlam)
  if(is.null(opts$alam)) opts$alam <- double(opts$nlam)
  if(is.null(opts$npass)) opts$npass <- integer(1)
  if(is.null(opts$jerr)) opts$jerr <- integer(1)
  
  M0 <- M
  U0 <- U
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2))
  n3 <- length(gam)
  nparams <- n1*n2*n3
  
  # The following lists save the corresponding objects for each tuning parameter, 
  B_l <- vector("list", nparams) # estimated B matrix
  beta_l <- vector("list", nparams) # estimated basis matrix
  step_l <- rep(NA_integer_, nparams) # iterations
  time_l <- rep(NA_real_, nparams) # running time
  rank_l <- rep(NA_integer_, nparams) # estimated rank
  s_l <- rep(NA_integer_, nparams) # estimated sparsity level
  
  # Count the number of converged matrices.
  nlam_cvg <- 0
  
  for(i in 1:n1){
    lambda1 <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){
        lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j,k])
        
        M <- M0 + gamma*diag(rep(1,ncol(M0)), ncol(M0),ncol(M0))
        
        # Initialize three matrices
        Bold <- matrix(0,dim(U0)[1], dim(U0)[2])
        Cold <- matrix(0,dim(U0)[1], dim(U0)[2])
        etaold <- matrix(0,dim(U0)[1], dim(U0)[2])
        
        # The MAIN loop of admm method
        step <- 0    
        start_time <- Sys.time()
        
        repeat{
          step <- step + 1
          
          # Update B
          U <- U0 - etaold + gamma * Cold
          out_B <- updateB(M, U, lambda1, opts)
          Bnew <- out_B$Bnew
          jerr <- out_B$jerr
          if(jerr != 0) break
          
          # Update C
          Cnew <- updateC(Bnew, lambda2, gamma, etaold)
          
          # Update eta (omega in SEAS algorithm)
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          # Code 1: success
          if(max(abs(Bnew - Cnew)) < eps){
            jerr <- 1
            break
          }
          # Code 404: then maximal iteration is reached
          if(step > maxit){
            jerr <- 404
            warning('Maximal iteration is reached.')
            break
          }
          Bold <- Bnew
          Cold <- Cnew
          etaold <- etanew
        }# End of repeat 
        end_time <- Sys.time()  # The time for each repeat
        time <- difftime(end_time, start_time, units = "secs")
        # Code < -10000: non-sparse matrix
        if(jerr < -10000){
          break
        }
        # Code 1: success, save the matrix and the related information.
        if(jerr==1){
          index <- (i-1)*n2*n3 + (j-1)*n2 + k
          nlam_cvg <- nlam_cvg + 1
          B_l[[index]] <- Bnew
          step_l[index] <- step
          time_l[index] <- time
          if(is.null(d)) rank <- rank_func(Cnew, thrd = eps)
          else rank <- d
          rank_l[index] <- rank
          # Cut and select the left singular vector of Bnew
          if(rank == 0){
            beta <- matrix(0, nrow(Bnew), ncol(Bnew))
          }else{
            tmp <- svd(Bnew)$u[,1:rank, drop = FALSE]
            vec <- as.vector(tmp)
            vec[abs(vec) < eps] <- 0
            beta <- matrix(vec, nrow(tmp), ncol(tmp))
          }
          beta_l[[index]] <- beta
          var_ind <- apply(beta, 1, function(x){any(x!=0)})
          s_l[index] <- sum(var_ind)
        }
      }# End of lambda2
      if(jerr < -10000) break
    }# End of gam
  }# End of lambda1
  return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg))
}

# Update B matrix in ADMM algorithm. Use the group-wise coordinate descent algorithm from 'msda' R package.
updateB <- function(M, U, lambda1, opts){
  U <- t(U)
  fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
  if(fit$jerr != 0){return(list(Bnew = NULL, jerr = fit$jerr))} # Code: non-zero, abnormal result.
  outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
  Bnew <- as.matrix(outlist$theta[[1]])
  list(Bnew = Bnew, jerr = fit$jerr)
}

# Update C matrix in ADMM algorithm.
updateC <- function(Bnew, lambda2, gamma, etaold){
  Btemp <- Bnew + 1/gamma * etaold
  svd_B <- svd(Btemp)
  lamtemp <- pmax(0, svd_B$d-lambda2/gamma)
  Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
  Cnew
}

# ------------------ revised functions from 'msda' package ---------------------- #
# Revise 'msda' function to accommodate other forms of M and U matrices.
# Some inputs are similar to the ones in 'seas' function. Please refer to 'msda' package documentation for more details of the arguments in 'msda' function.
# 
# Outputs:
# ========
# lambda: The tuning parameter sequence.
# theta: The list of estimated matrix.
# M: The M matrix from samples. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix from samples, which depends on the argument type.
# rank: The list of estimated rank for each matrix.
msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', FUN = NULL, lambda.factor=NULL, nlambda=100, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, verbose=FALSE, perturb=NULL){
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, type, FUN)
    M <- MU_out$M
    U <- MU_out$U
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax*2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) 
    diag(M) <- diag(M) + perturb
  H <- as.integer(dim(U)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam),alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, H)
  rank <- rep(NA_integer_, length(outlist$theta))
  for (i in 1:length(outlist$theta)){
    if(!is.null(outlist$theta[[i]])){
      rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
    }
  }
  if(is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  outlist <- list(lambda = outlist$lambda, theta = outlist$theta, M = M, U = U, rank = rank)
  class(outlist) <- c("msda")
  outlist
}

# Revise 'cv.msda' function to accommodate other forms of M and U matrices. We also add the optional argument 'fold' to pass the user-specified folds index.
# Some inputs are similar to the ones in 'cv.seas' function. Please refer to 'msda' package documentation for more details of arguments used in 'cv.msda' function.
# 
# Outputs:
# ========
# beta: The optimal estimated matrix.
# id: The index of the optimal tuning parameter.
# lambda: The lambda sequence.
# lam_max: The optimal tuning parameter.
# rank: The rank of the optimal estimated matrix.
# M: The M matrix based on the full data. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix based on the full data, which depends on the argument type.
# M_fold: The M matrix list based on each cross-validation data fold.
# U_fold: The U matrix list based on each cross-validation data fold.
cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', lambda.factor=NULL, nlambda=100, nfolds=5, foldid = NULL, lambda = NULL, FUN = NULL, maxit = 1e3, plot = FALSE){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  nobs <- nrow(x)
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  # Fit the model on the full data, obtain the lambda sequence.
  fit <- msda(x, y, yclass = yclass, type = type, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  
  # Cross-validation
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(foldid))
  }

  cv_out <- lapply(1:nfolds, function(k){
    x_train <- x[foldid!=k,,drop=FALSE]
    x_val <- x[foldid==k,,drop=FALSE]
    y_train <- y[foldid!=k]
    y_val <- y[foldid==k]
    yclass_train <- yclass[foldid!=k]
      
    fit_fold <- msda(x_train, y_train, yclass_train, type = type, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc(beta_fold, x_val, y_val)
    if(length(eval_fold) != length(lambda)){
      eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    }
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2)
  U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)){
   eval_all <- t(as.matrix(eval_all))
  }

  ## No matrix is converged in any fold
  if(all(is.na(eval_all))) return(NULL)
  
  cvm <- colMeans(eval_all, na.rm = TRUE)
  # The optimal lambda1
  id_max <- which.max(cvm)
  lam_max <- lambda[id_max]
  beta <- as.matrix(beta_l[[id_max]])
  
  # Recalculate the rank
  rank <- rank_func(beta, thrd = 1e-3)
  
  if(plot){ # If TRUE, plot the cv evaluation for each tuning parameter.
    dat <- data.frame(x = 1:length(cvm), y = cvm)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 1)+
      xlab("")+
      ylab("Distance correlation")+
      theme_bw()
    print(g)
  }
  
  list(beta = beta, id = id_max, lambda = lambda, lam_max = lam_max, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}
