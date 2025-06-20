data_gen_lmm_bsl <- function(n, m, p1, p2, q, k, mu, beta_t, beta1, beta2, mu_x,
                             Sigma_x, Omega_0, Sigma_d1, z, Sigma_e, theta){
  t <- rep(0:(m-1), n)
  
  X <- matrix(nrow=1, ncol=p1+p2)
  for (i in 1:n){
    X <- rbind(X,
               rep(1,m) %*% t(MASS::mvrnorm(n=1, mu_x, Sigma_x))
    )
  }
  X <- X[-1,]
  
  Xs <- list()
  for (r in 1:k){
    Delta <- matrix(nrow=1, ncol=p1+p2)
    for (i in 1:n){
      Delta <- rbind(Delta,
                     cbind(
                       rep(1,m) %*% t(MASS::mvrnorm(n=1, rep(0, p1), Sigma_d1)),
                       matrix(0, nrow=m, ncol=p2)
                     )
      )
    }
    Delta <- Delta[-1,]
    
    Xs[[r]] <- X + Delta
    Xs[[r]] <- as.matrix(Xs[[r]])
  }
  
  Sigma_d <- rbind( cbind(Sigma_d1, matrix(0, p1, p2)), matrix(0, p2, p1+p2) )
  
  b <- MASS::mvrnorm(n, rep(0, q), Omega_0)
  b <- c(t(b))
  
  z_ord <- order(z)
  z <- sort(z)
  Z <- list()
  for (i in 1:n){
    inds <- ((i-1)*m + 1) : (i*m)
    
    if (z[1] == -1) { # intercept
      if (length(z)>1 && z[2] == 0) { # intercept + time
        Z[[i]] <- cbind(
          rep(1,m),
          t[inds],
          X[inds, z[-(1:2)]]
        )
      } else { # intercept, no time
        Z[[i]] <- cbind(
          rep(1,m),
          X[inds, z[-1]]
        )
      }
    } else if (z[1] == 0){ # time, no intercept
      Z[[i]] <- cbind(
        t[inds],
        X[inds, z[-1]]
      )
    } else { # no intercept, no time
      Z[[i]] <- X[inds, z]
    }
    
    Z[[i]] <- Z[[i]][, z_ord]
  }
  Z <- bdiag(Z)
  
  e <- mvrnorm(n, rep(0, m), Sigma_e)
  e <- c(t(e))
  
  beta_x <- c(beta1, beta2)
  y <- c(as.matrix(
    mu + t*beta_t + X %*% beta_x + Z %*% b + e
  ))
  
  beta <- c(beta_t, beta_x)
  S <- which(beta != 0)
  
  #' The total number of unknown parameters is 2+p1+p2 (for intercept+time) + length(theta) + sig2e
  return (
    list(t=t, X=X, Xs=Xs, Delta=Delta, Z=Z, b=b, e=e, y=y,
         cluster=as.factor(rep(1:n, each=m)), beta_t=beta_t, beta_x=beta_x,
         beta1=beta1, beta2=beta2, beta=beta, S=S, mu=mu,
         mu_x=mu_x, Sigma_x=Sigma_x,
         Omega_0=Omega_0, Sigma_d1=Sigma_d1, Sigma_d=Sigma_d,
         Sigma_e=Sigma_e, theta=theta, n=n, m=m, N=n*m, p1=p1, p2=p2, p_x=p1+p2,
         p=1+p1+p2, q=q, k=k, d=length(theta), n_obs=n*m,
         n_params=2+p1+p2+1+length(theta))
  )
}

get_X_hat <- function(data, estimator){
  if (estimator=="_true"){
    Lambda_hat <- diag(data$p1+data$p2)
    X_hat <- data$X
  } else if (estimator=="_naive"){
    Lambda_hat <- diag(data$p1+data$p2)
    X_hat <- Reduce("+", data$Xs) / length(data$Xs)
  } else if (estimator=="_rck"){
    #' we can get the means without selecting only the baseline observations
    #' so long as there is an equal number of replicates (m) per subject. Then
    #' the mean including replicates is the same as mean with just baseline
    #' observations.
    Xs_bar <- Reduce("+", data$Xs) / length(data$Xs)
    Xs_bar_bsl <- Xs_bar[seq(from=1, by=data$m, length.out=data$n), ]
    mu_x_hat <- colMeans(Xs_bar_bsl)
    
    Sigma_xs <- data$Sigma_x + (data$Sigma_d / data$k)
    
    Lambda_hat <- solve(Sigma_xs) %*% data$Sigma_x
    X_hat <- rep(1,data$N) %*% t(mu_x_hat) %*% (diag(data$p_x) - Lambda_hat) +
      Xs_bar %*% Lambda_hat
  } else if (estimator=="_rcu"){
    #' we can get the means without selecting only the baseline observations
    #' so long as there is an equal number of replicates (m) per subject. Then
    #' the mean including replicates is the same as mean with just baseline
    #' observations.
    Xs_bar <- Reduce("+", data$Xs) / length(data$Xs)
    Xs_bar_bsl <- Xs_bar[seq(from=1, by=data$m, length.out=data$n), ]
    
    mu_x_hat <- colMeans(Xs_bar_bsl)
    Sigma_xs_bar_hat <- cov(Xs_bar_bsl)
    
    Sigma_delta_hat <- lapply(data$Xs, scale, center=TRUE, scale=FALSE)
    Sigma_delta_hat <- sapply(Sigma_delta_hat,
                              function(A, m, n) A[seq(from=1, by=m, length.out=n),],
                              m=data$m, n=data$n,
                              simplify="array")
    Sigma_delta_hat <- aperm(Sigma_delta_hat, c(3,2,1))
    Sigma_delta_hat <- apply(Sigma_delta_hat, MARGIN=3, FUN=cov,
                             simplify=FALSE)
    Sigma_delta_hat <- Reduce("+", Sigma_delta_hat) / data$n
    
    Sigma_x_hat <- (Sigma_xs_bar_hat - Sigma_delta_hat/data$k)
    
    Lambda_hat <- solve(Sigma_xs_bar_hat) %*% Sigma_x_hat
    X_hat <- rep(1,data$N) %*% t(mu_x_hat) %*% (diag(data$p_x) - Lambda_hat) +
      Xs_bar %*% Lambda_hat
  } else {
    stop("invalid method '", estimator, "'")
  }
  
  return ( X_hat )
}

get_OmegaBias <- function(data, beta, estimator){
  beta_t <- beta[1]
  beta_x <- beta[2:(1+data$p_x)]
  
  if ( estimator == "_rck" ){
    Xs_bar <- Reduce("+", data$Xs) / length(data$Xs)
    # valid if same number of replicates per subject
    Xs_bar_bsl <- Xs_bar[seq(from=1, by=data$m, length.out=data$n), ]
    mu_x_hat <- colMeans(Xs_bar_bsl)
    
    Sigma_xs <- data$Sigma_x + (data$Sigma_d / data$k)
    
    Lambda_hat <- solve(Sigma_xs) %*% data$Sigma_x
    
    Omega_tau_0 <- matrix(c(
      t(beta_x) %*% data$Sigma_x %*%
        (diag(data$p_x) - Lambda_hat) %*% beta_x
    ),
    nrow=data$q, ncol=data$q)
  } else if ( estimator == "_rcu" ){
    Xs_bar <- Reduce("+", data$Xs) / length(data$Xs)
    # valid if same number of replicates per subject
    Xs_bar_bsl <- Xs_bar[seq(from=1, by=data$m, length.out=data$n), ]
    mu_x_hat <- colMeans(Xs_bar_bsl)
    Sigma_xs_bar_hat <- cov(Xs_bar_bsl)
    
    Sigma_delta_hat <- lapply(data$Xs, scale, center=TRUE, scale=FALSE)
    Sigma_delta_hat <- sapply(Sigma_delta_hat,
                              function(A, m, n) A[seq(from=1, by=m, length.out=n),],
                              m=data$m, n=data$n,
                              simplify="array")
    Sigma_delta_hat <- aperm(Sigma_delta_hat, c(3,2,1))
    Sigma_delta_hat <- apply(Sigma_delta_hat, MARGIN=3, FUN=cov,
                             simplify=FALSE)
    Sigma_delta_hat <- Reduce("+", Sigma_delta_hat) / data$n
    
    Sigma_x_hat <- (Sigma_xs_bar_hat - Sigma_delta_hat/data$k)
    
    Lambda_hat <- solve(Sigma_xs_bar_hat) %*% Sigma_x_hat
    
    Omega_tau_0 <- matrix(c(
      t(beta_x) %*% Sigma_x_hat %*%
        (diag(data$p_x) - Lambda_hat) %*% beta_x
    ),
    nrow=data$q, ncol=data$q)
  } else {
    Omega_tau_0 <- 0
  }
  
  return ( Omega_tau_0 )
}