icaimax <-
  function(X, nc, center = TRUE, maxit = 100, tol = 1e-6, Rmat = diag(nc), 
           alg = "newton", fun = "tanh", signs = rep(1, nc), signswitch = TRUE, 
           rate = 1, rateanneal = NULL){
    ###### ICA via (Fast and Robust) Information-Maximization
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 3, 2022
    
    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc < 1) stop("Must set nc >= 1 component") 
    maxit <- as.integer(maxit[1])
    if(maxit < 1) stop("Must set maxit >= 1 iteration")
    tol <- tol[1]
    if(tol <= 0) stop("Must set tol > 0") 
    if(nc > min(nobs, nvar)) stop("Too many components. Set nc <= min(dim(X))") 
    if(nrow(Rmat) != nc | ncol(Rmat) != nc) stop("Input 'Rmat' must be nc-by-nc rotation matrix.")
    fun <- fun[1]
    if(fun == "ext"){
      signs <- sign(signs)
      if(length(signs) != nc){ stop("Input 'signs' must be have length equal to 'nc' input.") }
    } else {
      signs <- NA
      signswitch <- FALSE
    }
    alg <- alg[1]
    if(alg == "gradient"){
      rate <- rate[1]
      if(rate <= 0) stop("Must set 'rate' greater than 0")
      if(!is.null(rateanneal[1])){
        if(length(rateanneal) != 2L) stop("Input 'rateanneal' should be two-element vector")
        if(rateanneal[1] <= 0 || rateanneal[1] >= 90) stop("Input 'rateanneal[1]' should be in range (0, 90)")
        if(rateanneal[2] <= 0 || rateanneal[2] > 1){ stop("Input 'rateanneal[2]' should be in range (0, 1]") }
        ralog <- TRUE
      } else {
        ralog <- FALSE
      }
    }
    
    ### center and whiten
    if(center) X <- scale(X, scale = FALSE)
    if(nobs >= nvar){
      xeig <- eigen(crossprod(X) / nobs, symmetric = TRUE)
    } else {
      xeig <- eigen(tcrossprod(X) / nobs, symmetric = TRUE)
    } # end if(nobs >= nvar)
    nze <- sum(xeig$values > xeig$values[1] * .Machine$double.eps)
    if(nze < nc){
      warning("Numerical rank of X is less than requested number of components (nc).\nNumber of components has been redefined as rank(X) = ",nc)
      nc <- nze
      Rmat <- diag(nc)
    }
    Dmat <- sdiag(sqrt(xeig$values[1:nc]))
    if(nobs >= nvar){
      Mprt <- tcrossprod(Dmat, xeig$vectors[, 1:nc, drop = FALSE])
      diag(Dmat) <- 1 / diag(Dmat)
      Pmat <- xeig$vectors[, 1:nc, drop = FALSE] %*% Dmat
      Xw <- X %*% Pmat   # whitened data
    } else {
      Mprt <- crossprod(xeig$vectors[, 1:nc, drop = FALSE], X) / sqrt(nobs)
      diag(Dmat) <- 1 / diag(Dmat)^2
      Pmat <- crossprod(Mprt, Dmat)
      Xw <- xeig$vectors[, 1:nc, drop = FALSE] * sqrt(nobs)   # whitened data
    } # end if(nobs >= nvar)
    
    ### check if nc=1
    if(nc == 1L){
      res <- list(S = Xw, M = Mprt, W = t(Pmat), Y = Xw, Q = t(Pmat), 
                  R = matrix(1), vafs = nobs * sum(Mprt^2) / sum(X^2), 
                  iter = NA, alg = alg[1], fun = fun[1], signs = signs,
                  rate = rate, converged = TRUE)
      class(res) <- "icaimax"
      return(res)
    }
    
    ### which nonlinearity
    if(fun == "log"){
      fun1d <- function(x, sgn = 1){ 2 / (1 + exp(-x)) - 1 }
      fun2d <- function(x, sgn = 1){ 1 / (cosh(x) + 1) }
    } else if(fun == "ext"){
      fun1d <- function(x, sgn = 1){ x + tanh(x) %*% sdiag(sgn) }
      fun2d <- function(x, sgn = 1){ 1 + (1 - tanh(x)^2) %*% sdiag(sgn) }
    } else {
      fun1d <- function(x, sgn = 1){ tanh(x) }
      fun2d <- function(x, sgn = 1){ 1 - tanh(x)^2 }
    }
    
    ### which algorithm
    if(alg[1] == "gradient"){
      
      # gradient descent
      iter <- 0
      vtol <- 1
      while(vtol > tol && iter < maxit){
        # update all components
        smat <- Xw %*% Rmat
        if(signswitch) signs <- sign(colMeans((cosh(smat)^-2) - tanh(smat) * smat)) 
        rnew <- Rmat - rate * crossprod(Xw / nobs, fun1d(smat, signs))
        # orthgonalize
        rsvd <- svd(rnew)
        rnew <- tcrossprod(rsvd$u, rsvd$v)
        # check for convergence
        vtol <- 1 - min(abs(colSums(Rmat * rnew)))
        iter <- iter + 1
        Rmat <- rnew
        if(ralog && ((acos(1 - vtol) * 180 / pi) < rateanneal[1])) rate <- rate * rateanneal[2]
      } # end while(vtol>tol && iter<maxit)
      
    } else {
      
      # Newton iteration
      iter <- 0
      vtol <- 1
      while(vtol > tol && iter < maxit){
        # update all components
        smat <- Xw %*% Rmat
        if(signswitch) signs <- sign(colMeans((cosh(smat)^-2) - tanh(smat) * smat)) 
        Hmat <- matrix(colMeans(fun2d(smat, signs)), nrow = nc, ncol = nc, byrow = TRUE)
        rnew <- Rmat - crossprod(Xw / nobs, fun1d(smat, signs)) / Hmat
        # orthgonalize
        rsvd <- svd(rnew)
        rnew <- tcrossprod(rsvd$u, rsvd$v)
        # check for convergence
        vtol <- 1 - min(abs(colSums(Rmat * rnew)))
        iter <- iter + 1
        Rmat <- rnew
      } # end while(vtol>tol && iter<maxit)
      
    } # end if(alg=="gradient")
    
    ### sort according to vafs
    M <- crossprod(Rmat, Mprt)
    vafs <- rowSums(M^2)
    ix <- sort(vafs, decreasing = TRUE, index.return = TRUE)$ix
    M <- M[ix,]
    Rmat <- Rmat[,ix]
    vafs <- nobs * vafs[ix]  / sum(X^2)
    
    ### return results
    res <- list(S = Xw %*% Rmat, M = t(M), W = t(Pmat %*% Rmat), Y = Xw,
                Q = t(Pmat), R = Rmat, vafs = vafs, iter = iter,
                alg = alg[1], fun = fun[1], signs = signs, rate = rate,
                converged = ifelse(vtol <= tol, TRUE, FALSE))
    class(res) <- "icaimax"
    return(res)
    
  }

print.icaimax <- 
  function(x, ...){
    nc <- length(x$vafs)
    cat("\nInfomax ICA with", nc, ifelse(nc == 1L, "component", "components"), "\n\n")
    cat("   converged: ", x$converged,"  (", x$iter, " iterations) \n", sep = "")
    cat("   r-squared:", sum(x$vafs), "\n")
    cat("   algorithm:", x$alg, "\n")
    cat("    function:", x$fun, "\n\n")
  }