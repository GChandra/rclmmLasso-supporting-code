#' Code for LMM Lasso algorithm

# soft-thresholding algorithm
st <- function(x, alpha){
  st_x <- rep(0, length(x))
  
  for (i in 1:length(x)){
    if (x[i] > alpha){ st_x[i] <- x[i] - alpha }
    else if (x[i] < -alpha){ st_x[i] <- x[i] + alpha }
    else { st_x[i] <- 0 }
  }
  
  return ( st_x )
}

#' Deviance function (twice negative log-likelihood) parametrized by variance
#' components.
devtheta <- function(theta, y, X, Z, L, Lambda, sig2e, uk, beta, n, N, p, q,
                     R_XtR_X=NULL, REML=TRUE){
  u <- uk
  if (REML & is.null(R_XtR_X)){
    REML <- FALSE
    message("R_X matrix not provided; cannot do REML optimization.
            Switching to normal optimization.")
  }
  
  Lambda@x[] <- mapping(theta, n, q, init=FALSE)
  L <- update_L(L, Lambda, Z, init=FALSE)
  
  r2theta <- r2theta(y, X, Z, Lambda, u, beta)
  
  if (REML){
    return(
      2*determinant(L, sqrt=TRUE, logarithm=TRUE)$modulus +
        determinant(R_XtR_X, logarithm=TRUE)$modulus +
        (N-p)*(1 + log(2*pi*r2theta/(N-p)))
    )
  } else {
    return(
      2*determinant(L, sqrt=TRUE, logarithm=TRUE)$modulus + N*(1 + log(2*pi*r2theta/N))
    )
  }
}

#' Residual function
r2theta <- function(y, X, Z, Lambda, u, beta){
  return (
    norm(y - X%*%beta - Z%*%Lambda%*%u, "2")^2 + norm(u, "2")^2
  )
}

#' Mapping of variance components to random effects covariance matrix
mapping <- function(theta, n, q, init=TRUE){
  if ( init ){
    Lambda_0 <- matrix(0, q, q)
    Lambda_0[lower.tri(Lambda_0, diag=TRUE)] <- theta
    Lambda <- bdiag( rep(list(Lambda_0), n) )
    
    return( Lambda )
  } else {
    return ( rep(theta, n) )
  }
}

update_L <- function(L, Lambda, Z, init=TRUE){
  if ( init ){
    LLt <- t(Z %*% Lambda) %*% (Z %*% Lambda) + diag(dim(Z)[2])
    return( Cholesky(LLt, LDL=FALSE) )
  } else {
    return( Matrix::update(L, t(Z%*%Lambda), mult=1) )
  }
}

#' LMM fitting algorithm
lmm_lasso <- function(y, X, Z,
                      theta_0=NULL, beta_0=NULL, sig2e_0=NULL,
                      p, q, d, n, N,
                      REML=TRUE, BIC=TRUE, lambda, penalty.factor=NULL,
                      eta=1e-4, eta_gs=1e-3, maxits=1000,
                      standardize=TRUE, intercept=TRUE){
  # center data
  y0 <- as.matrix(y) - mean(y)*intercept
  X0 <- scale(X, center=intercept, scale=standardize)
  
  if (is.null(penalty.factor)){
    penalty.factor <- rep(1, p)
  } else {
    # rescale factors to add to p, as in glmnet
    penalty.factor <- penalty.factor * p / sum(penalty.factor)
  }
  
  if (is.null(theta_0) | is.null(beta_0) | is.null(sig2e_0)){
    beta_k <- rep(0, p)
    theta_k <- rep(1, d)
    sig2e_k <- 1
  } else {
    beta_k <- beta_0
    theta_k <- theta_0
    sig2e_k <- sig2e_0
  }
  
  # create lower bound vector for nloptr (lower-triangular reference is based on
  # the choice in the mapping() function).
  lb <- matrix(-Inf, nrow=q, ncol=q)
  diag(lb) <- 0
  lb <- lb[lower.tri(lb, diag=TRUE)]
  
  # create covariance Cholesky factor
  Lambda_k <- mapping(theta_k, n, q, init=TRUE)
  L <- update_L(L=NULL, Lambda=Lambda_k, Z=Z, init=TRUE)
  
  u_k <- rep(0, n*q)
  d_k <- c(u_k, beta_k)
  
  d_old <- d_k + 2*eta
  
  k <- 0
  while ( (max(abs(d_k-d_old)) > eta) & (k < maxits) ) {
    d_old <- d_k
    
    # update covariance Cholesky factor
    Lambda_k@x[] <- mapping(theta_k, n, q, init=FALSE)
    L <- update_L(L=L, Lambda=Lambda_k, Z=Z, init=FALSE)
    
    # Solve normal equations
    R_ZX <- solve( L, crossprod(Z%*%Lambda_k, X0), system="L" )
    R_XtR_X <- as( crossprod(X0) - crossprod(R_ZX), "dpoMatrix" )
    
    # outer loop for convergence
    # set active-set
    J <- which(beta_k != 0)
    if ( (k < 20) || (k %% 200 == 0) ){
      J <- 1:p
    }
    
    d_k_i <- c(u_k, beta_k)
    d_old_i <- d_k_i + 2*eta_gs
    while(max(abs(d_k_i-d_old_i)) > eta_gs){
      d_old_i <- d_k_i
      
      u_k <- c(as.matrix(
        solve(L, crossprod(Z%*%Lambda_k, y0) - crossprod(Z%*%Lambda_k, X0%*%beta_k),
              system="LDLt")
      ))
      
      max_diff <- 2*eta_gs
      while (max_diff > eta_gs){
        max_diff <- -Inf
        
        for (j in J){
          beta_j_new <- st( t(X0[,j]) %*%
                              (y0 - X0[,-j]%*%beta_k[-j] - Z%*%Lambda_k%*%u_k) /
                              norm(X0[,j], "2")^2,
                            lambda * penalty.factor[j] / norm(X0[,j], "2")^2)
          
          if (abs(beta_j_new - beta_k[j]) > max_diff){
            max_diff <- abs(beta_j_new - beta_k[j])
          }
          
          beta_k[j] <- beta_j_new
        }
      }
      
      d_k_i <- c(u_k, beta_k)
    }
    
    d_k <- d_k_i
    
    # Update residual
    if (REML){
      sig2e_k <- r2theta(y0, X0, Z, Lambda_k, u_k, beta_k) / (N-p)
    } else {
      sig2e_k <- r2theta(y0, X0, Z, Lambda_k, u_k, beta_k) / N
    }
    
    # optimization step
    opt_dev <- nloptr(x0=theta_k, eval_f=devtheta, eval_grad_f=NULL,
                      lb=lb,
                      opts=list("algorithm"="NLOPT_LN_BOBYQA",
                                "xtol_rel"=1.0e-4),
                      y=y0, X=X0, Z=Z, L=L, Lambda=Lambda_k,
                      sig2e=sig2e_k, uk=u_k,
                      beta=beta_k,
                      n=n, N=N, p=p, q=q,
                      R_XtR_X=R_XtR_X, REML=REML)
    # adding 1e-10 to avoid issues with zeroing out elements I need to access.
    theta_k <- opt_dev$solution
    theta_k <- theta_k + 1e-10*sign(theta_k)
    
    k <- k + 1
  }
  
  if (k == maxits){ warning("algorithm failed to converge, reached max iterations.") }
  
  # update covariance Cholesky factor and u_k once more
  Lambda_k@x[] <- mapping(theta_k, n, q, init=FALSE)
  L <- update_L(L=L, Lambda=Lambda_k, Z=Z, init=FALSE)
  # Solve normal equations
  R_ZX <- solve( L, crossprod(Z%*%Lambda_k, X0), system="L" )
  R_XtR_X <- as( crossprod(X0) - crossprod(R_ZX), "dpoMatrix" )
  u_k <- c(as.matrix(
    solve(L, crossprod(Z%*%Lambda_k, y0) - crossprod(Z%*%Lambda_k, X0%*%beta_k),
          system="LDLt")
  ))
  
  if (standardize){
    beta_k <- beta_k / attr(X0, "scaled:scale")
  }
  
  if (BIC) {
    ic <- opt_dev$objective + log(N)*(sum(beta_k != 0) + d)
  } else {
    ic <- opt_dev$objective + 2*(sum(beta_k != 0) + d)
  }
  
  mu_k <- 0
  if (intercept){
    mu_k <- mean(y - X %*% beta_k - Z %*% Lambda_k %*% u_k)
  }
  Omega_0_opt <- as.matrix(Lambda_k[1:q,1:q]) %*% t(as.matrix(Lambda_k[1:q,1:q])) * sig2e_k
  
  return (list(
    eta_opt=max(abs(d_k-d_old)),
    ic_opt=ic,
    theta_opt=theta_k,
    Omega_0_opt=Omega_0_opt,
    mu_opt=mu_k,
    beta_opt=beta_k,
    u_opt=u_k,
    sig2e_opt=sig2e_k,
    iter=k
  ))
}

# ... - extra parameters for lmer function
ic.lmm_lasso <- function(formula, data,
                          p, q, d, n, N,
                          theta_0=NULL, beta_0=NULL, sig2e_0=NULL,
                          lambda.max=NULL, nlambda=20, lambda.min.ratio=0.001,
                          penalty.factor=NULL, REML=TRUE, BIC=TRUE,
                          eta=1e-4, eta_gs=5e-3, maxits=5e4,
                          standardize=TRUE, progress=FALSE,
                          getInit=FALSE, ...){
  if (is.null(penalty.factor)){
    penalty.factor <- rep(1, p)
  } else {
    # rescale factors to sum to p, as in glmnet
    penalty.factor <- penalty.factor * p / sum(penalty.factor)
  }
  
  formula <- formula(formula)
  
  frm <- lFormula(eval(formula), data)
  Z <- t(frm$reTrms$Zt)
  
  intercept <- all(frm$X[,1]==1)
  y <- frm$fr[,1]
  if (intercept){
    X <- frm$X[,-1]
  } else {
    X <- frm$X
  }
  
  y0 <- y - mean(y)*intercept
  X0 <- scale(X, center=intercept, scale=standardize)
  
  if (is.null(lambda.max)){
    if(any(penalty.factor==0)){
      fit <- lmer(
        paste0(frm$formula[[2]], " ~ ", as.numeric(intercept), " + ",
               paste(colnames(frm$X)[-1][penalty.factor==0], collapse=" + "), " + ",
               paste("(", findbars(formula), ")", collapse=" + ")
        ),
        data,
        ...
      )
    } else {
      fit <- lmer(
        paste0(frm$formula[[2]], " ~ ", as.numeric(intercept), " + ",
               paste("(", findbars(formula), ")", collapse=" + ")
        ),
        data,
        ...
      )
    }
    
    beta_0 <- rep(0, p)
    beta_0[penalty.factor==0] <- fixef(fit)[-1]
    sig2e_0 <- as.data.frame(VarCorr(fit))$vcov[d+1]
    theta_0 <- c(vech(t(chol(as.matrix(Matrix::bdiag(VarCorr(fit))))))) / sqrt(sig2e_0)
    
    lambda.max.ub <- max(abs( t(X0) %*% y0 )) * sum(penalty.factor) / p
    fit <- lmm_lasso(y, X, Z,
                     theta_0, beta_0, sig2e_0,
                     p, q, d, n, N,
                     REML=REML, BIC=BIC,
                     lambda=lambda.max.ub, penalty.factor=penalty.factor,
                     eta=eta, eta_gs=eta_gs, maxits=maxits,
                     standardize=standardize, intercept=intercept)
    theta_0 <- fit$theta_opt
    beta_0 <- fit$beta_opt
    sig2e_0 <- fit$sig2e_opt
    
    # create covariance Cholesky factor
    Lambda_0 <- mapping(theta_0, n, q)
    
    lambda.max <- max(abs( t(X0) %*% (y0 -
                                        X0%*%beta_0 - Z%*%Lambda_0%*%fit$u_opt) )) *
      sum(penalty.factor) / p
    
    if (getInit){
      return (list(
        lambda.max=lambda.max,
        theta_0=theta_0,
        beta_0=beta_0,
        sig2e_0=sig2e_0
      ))
    }
  }
  
  lambda.seq <- exp(
    seq( log(lambda.max), log(lambda.min.ratio*lambda.max), length.out=nlambda )
  )
  
  if (progress){
    start <- proc.time()
  }
  S_hat <- matrix(NA, nrow=nlambda, ncol=p)
  ic_min <- Inf
  for (lambda in lambda.seq){
    if (progress){
      print(lambda)
    }
    
    fit <- lmm_lasso(y, X, Z,
                     theta_0, beta_0, sig2e_0,
                     p, q, d, n, N,
                     REML=REML, BIC=BIC,
                     lambda=lambda, penalty.factor=penalty.factor,
                     eta=eta, eta_gs=eta_gs, maxits=maxits,
                     standardize=standardize, intercept=intercept)
    theta_0 <- fit$theta_opt
    beta_0 <- fit$beta_opt
    
    if (progress==2){
      beta_p <- beta_0
      names(beta_p) <- NULL
      print(beta_p)
      print(theta_0)
      print(fit$ic_opt)
    }
    if (standardize){
      beta_0 <- beta_0 * attr(X0, "scaled:scale")
    }
    
    sig2e_0 <- fit$sig2e_opt
    
    S_hat[ which(lambda.seq==lambda) , ] <- as.numeric(beta_0 != 0)
    
    if (fit$ic_opt < ic_min){
      ic_min <- fit$ic_opt
      
      theta_opt <- fit$theta_opt
      mu_opt <- fit$mu_opt
      beta_opt <- fit$beta_opt
      u_opt <- fit$u_opt
      sig2e_opt <- fit$sig2e_opt
      Omega_0_opt <- fit$Omega_0_opt
      lambda_opt <- lambda
    }
  }
  if (progress){
    end <- proc.time()
    print(paste( round((end-start)[3], 2), "seconds" ))
  }
  
  return (list(
    ic_min=ic_min,
    mu_opt=mu_opt,
    beta_opt=beta_opt,
    u_opt=u_opt,
    S_hat=S_hat,
    sig2e_opt=sig2e_opt,
    theta_opt=theta_opt,
    Omega_0_opt=Omega_0_opt,
    lambda_opt=lambda_opt
  ))
}
