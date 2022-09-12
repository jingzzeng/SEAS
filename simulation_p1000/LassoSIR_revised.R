# Revise the 'LassoSIR' function from R package 'LassoSIR':
# (1) Record the computation time. 
# (2) Add user-specified fold id ("foldid").
# (3) Add user-specified lambda sequence ("lam").
LassoSIR_revised <- function (X, Y, H = 0, choosing.d = "automatic", solution.path = FALSE, categorical = FALSE, nfolds = 10, foldid = NULL, screening = TRUE, no.dim = 0, lam = NULL) 
{
  start_time <- Sys.time()
  if (no.dim != 0) 
    choosing.d = "given"
  if ((categorical == FALSE) & (H == 0)) {
    H <- (function() {
      H <- readline("For the continuous response, please choose the number of slices:   ")
      H <- as.numeric(unlist(strsplit(H, ",")))
      return(dim)
    })()
  }
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (categorical == FALSE) {
    ORD <- order(Y)
    X <- X[ORD, ]
    Y <- Y[ORD]
    ms <- array(0, n)
    m <- floorH)
    c <- n%%H
    M <- matrix(0, nrow = H, ncol = n)
    if (c == 0) {
      M <- diag(H) %x% matrix(1, nrow = 1, ncol = m
      ms <- m + ms
    }
    else {
      for (i in 1:c) {
        M[i, ((m + 1) * (i - 1) + 1):((m + 1) * i)] <-(m + 
                                                            1)
        ms[((m + 1) * (i - 1) + 1):((m + 1) * i)] <- m
      }
      for (i in (c + 1):H) {
        M[i, ((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                     1) * c + (i - c) * m)] <-m
        ms[((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                   1) * c + (i - c) * m)] <- m - 1
      }
    }
    if (screening == TRUE) {
      x.sliced.mean <- M %*% X
      sliced.variance <- apply(x.sliced.mean, 2, var)
      keep.ind <- sort(order(sliced.variance, decreasing = TRUE)[1:n])
    }
    else {
      keep.ind <- c(1:p)
    }
    X <- X[, keep.ind]
    X.H <- matrix(0, nrow = H, ncol = dim(X)[2])
    grand.mean <- matrix(apply(X, 2, mean), nrow = 1, ncol = dim(X)[2])
    X.stand.ord <- X - grand.mean %x% matrix(1, nrow = dim(X)[1], 
                                             ncol = 1)
    X.H <- M %*% X.stand.ord
  }
  else {
    ms <- array(0, n)
    Y.unique <- unique(Y)
    H <- length(Y.unique)
    ORD <- which(Y == Y.unique[1])
    nH <- sum(Y == Y.unique[1])
    ms[1:nH] <- nH
    for (i in 2:H) {
      ORD <- c(ORD, which(Y == Y.unique[i]))
      nH <- c(nH, sum(Y == Y.unique[i]))
      ms[(sum(nH[1:(i - 1)]) + 1):sum(nH[1:i])] <- nH[i]
    }
    X <- X[ORD, ]
    M <- matrix(0, nrow = H, ncol = n)
    M[1, 1:nH[1]] <-nH[1]
    for (i in 2:H) M[i, (sum(nH[1:(i - 1)]) + 1):sum(nH[1:i])] <-nH[i]
    if (screening == TRUE) {
      x.sliced.mean <- M %*% X
      sliced.variance <- apply(x.sliced.mean, 2, var)
      keep.ind <- sort(order(sliced.variance, decreasing = TRUE)[1:n])
    }
    else {
      keep.ind <- c(1:p)
    }
    X <- X[, keep.ind]
    X.H <- matrix(0, nrow = H, ncol = dim(X)[2])
    grand.mean <- matrix(apply(X, 2, mean), nrow = 1, ncol = dim(X)[2])
    X.stand.ord <- X - grand.mean %x% matrix(1, nrow = dim(X)[1], 
                                             ncol = 1)
    X.H <- M %*% X.stand.ord
  }
  svd.XH <- svd(X.H, nv = p)
  res.eigen.value <- array(0, p)
  res.eigen.value[1:dim(X.H)[1]] <- (svd.XH$d)H
  if (choosing.d == "manual") {
    plot(c(1:p), res.eigen.value, ylab = "eigen values")
    no.dim <- (function() {
      dim <- readline("Choose the number of directions:   ")
      dim <- as.numeric(unlist(strsplit(dim, ",")))
      return(dim)
    })()
  }
  # Record time
  end_time <- Sys.time()
  time1 <- difftime(end_time, start_time, units = "secs")
  
  if (choosing.d == "automatic") {
    beta.hat <- array(0, c(p, min(p, H)))
    Y.tilde <- array(0, c(n, min(p, H)))
    for (ii in 1:min(p, H)) {
      eii <- matrix(0, nrow = dim(svd.XH$v)[2], ncol = 1)
      eii[ii] <- 1
      eigen.vec <- solve(t(svd.XH$v), eii)
      Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.v(res.eigen.value[ii]) * 
        matrixms, nrow = n, ncol = 1)
    }
    mus <- array(0, min(p, H))
    for (ii in 1:min(p, H)) {
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[, ii], nfolds = nfolds, foldid = foldid, lambda = lam)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      mus[ii] <- lambda
      lars.fit <- glmnet(X.stand.ord, Y.tilde[, ii], lambda = lambda)
      beta.hat[keep.ind, ii] <- as.double(lars.fit$beta)
    }
    temp.2 <- sqrt(apply(beta.hat^2, 2, sum)) * res.eigen.value[1:H]
    temp <- temptemp.2[1]
    res.kmeans <- kmeans(temp, centers = 2)
    no.dim <- min(sum(res.kmeans$cluster == 1), sum(res.kmeans$cluster == 
                                                      2))
  }
  start_time <- Sys.time()
  beta.hat <- array(0, c(p, no.dim))
  Y.tilde <- array(0, c(n, no.dim))
  for (ii in 1:no.dim) {
    eii <- matrix(0, nrow = dim(t(svd.XH$v))[2], ncol = 1)
    eii[ii] <- 1
    eigen.vec <- solve(t(svd.XH$v), eii)
    Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.v(res.eigen.value[ii]) * 
      matrixms, nrow = n, ncol = 1)
  }
  # Record time
  end_time <- Sys.time()
  time2 <- difftime(end_time, start_time, units = "secs")
  
  if (solution.path == FALSE) {
    mus <- array(0, no.dim)
    if (no.dim == 1) {
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde, nfolds = nfolds, foldid = foldid, lambda = lam)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      start_time <- Sys.time()
      lars.fit <- glmnet(X.stand.ord, Y.tilde, lambda = lambda)
      beta.hat[keep.ind] <- as.double(lars.fit$beta)
      # Record time
      end_time <- Sys.time()
      time3 <- difftime(end_time, start_time, units = "secs")
    }
    else {
      time3 <- c()
      for (ii in 1:no.dim) {
        lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[, ii], nfolds = nfolds, foldid = foldid, lambda = lam)
        ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
        if (ind == 1) 
          ind <- 2
        lambda <- lars.fit.cv$lambda[ind]
        mus[ii] <- lambda
        start_time <- Sys.time()
        lars.fit <- glmnet(X.stand.ord, Y.tilde[, ii], lambda = lambda)
        beta.hat[keep.ind, ii] <- as.double(lars.fit$beta)
        # Record time
        end_time <- Sys.time()
        time3 <- c(time3, difftime(end_time, start_time, units = "secs"))
      }
    }
    time <- time1 + time2 + sum(time3) # We do not include the time for dimension selection and tuning parameter selection.
    list(beta = beta.hat, eigen.value = res.eigen.value, 
         no.dim = no.dim, H = H, categorical = categorical, time = time,           lam = lars.fit.cv$lambda)
  }
  else {
    lars.fit.all <- list()
    for (ii in 1:no.dim) {
      lars.fit.all[[ii]] <- glmnet(X.stand.ord, Y.tilde[, ii])
    }
    lars.fit.all
  }
}