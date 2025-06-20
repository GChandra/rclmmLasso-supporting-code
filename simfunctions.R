run_sim <- function(data, estimator, runs, of, nlambda=20, REML=TRUE, eta=1e-4,
                    eta_gs=5e-3, maxits=5e4, standardize=TRUE, progress=FALSE,
                    parallel=FALSE, mc.cores=NULL){
  if(parallel){
    cat("running condition in parallel\n", file=of, append=TRUE)
    return(
      run_sim_p(data=data, estimator=estimator, runs=runs, of=of, nlambda=nlambda,
                REML=REML, eta=eta, eta_gs=eta_gs, maxits=maxits,
                standardize=standardize, progress=progress, mc.cores=mc.cores)
    )
  } else {
    cat("running condition serially\n", file=of, append=TRUE)
    return(
      run_sim_s(data=data, estimator=estimator, runs=runs, of=of, nlambda=nlambda,
                REML=REML, eta=eta, eta_gs=eta_gs, maxits=maxits,
                standardize=standardize, progress=progress)
    )
  }
}

#' run simulations in parallel across runs
run_sim_p <- function(data, estimator, runs, of, nlambda, REML, eta, eta_gs,
                      maxits, standardize, progress, mc.cores){
  runs_start <- proc.time()
  
  res_t <- mclapply(data, run_iter, estimator=estimator, nlambda=nlambda,
                  REML=REML, eta=eta, eta_gs=eta_gs, maxits=maxits,
                  standardize=standardize, progress=progress, mc.cores=mc.cores)
  res_t <- simplify2array(res_t)
  
  runs_time <- proc.time()
  cat("runs completed in: ", round((runs_time-runs_start)[3], 2),
      " seconds\n", file=of, append=TRUE)
  
  l2_mu <- res_t[1,]
  l2 <- res_t[2,]
  l2S <- res_t[3,]
  lF_Omega <- res_t[4,]
  fpr_opt <- res_t[5,]
  tpr_opt <- res_t[6,]
  fpr <- res_t[7:(7+nlambda-1),]
  tpr <- res_t[(7+nlambda):(7+2*nlambda-1),]
  
  res <- c(l2_mu=mean(l2_mu), l2_mu_sd=sd(l2_mu),
           l2=mean(l2), l2_sd=sd(l2),
           l2S=mean(l2S), l2S_sd=sd(l2S),
           lF_Omega=mean(lF_Omega), lF_Omega_sd=sd(lF_Omega),
           fpr_opt=mean(fpr_opt), fpr_opt_sd=sd(fpr_opt),
           tpr_opt=mean(tpr_opt), tpr_opt_sd=sd(tpr_opt),
           fpr=rowMeans(fpr), tpr=rowMeans(tpr))
  res <- round(res, 3)
  
  return (res)
}

#' run simulations serially (not in parallel)
run_sim_s <- function(data, estimator, runs, of, nlambda, REML, eta, eta_gs,
                      maxits, standardize, progress){
  fpr <- tpr <- matrix(NA, nrow=runs, ncol=nlambda)
  l2_mu <- l2 <- l2S <- lF_Omega <- fpr_opt <- tpr_opt <- rep(0, runs)
  
  runs_start <- proc.time()
  upd_start <- runs_start
  
  for (r in 1:runs){
    cat("run:", r, "\n", file=of, append=TRUE)
    if (r %% 25 == 0){
      upd_time <- proc.time()
      cat("progress update: run ", r, ": ",
          round((upd_time-upd_start)[3], 2), " seconds since last update\n",
          file=of, append=TRUE)
      upd_start <- proc.time()
    }
    
    res <- run_iter(data=data[[r]], estimator=estimator, nlambda=nlambda,
                    REML=REML, eta=eta, eta_gs=eta_gs, maxits=maxits,
                    standardize=standardize, progress=progress)
    l2_mu[r] <- res[1]
    l2[r] <- res[2]
    l2S[r] <- res[3]
    lF_Omega[r] <- res[4]
    fpr_opt[r] <- res[5]
    tpr_opt[r] <- res[6]
    fpr[r,] <- res[7:(7+nlambda-1)]
    tpr[r,] <- res[(7+nlambda):(7+2*nlambda-1)]
  }
  
  runs_time <- proc.time()
  cat("runs completed in: ", round((runs_time-runs_start)[3], 2),
      " seconds\n", file=of, append=TRUE)
  
  res <- c(l2_mu=mean(l2_mu), l2_mu_sd=sd(l2_mu),
           l2=mean(l2), l2_sd=sd(l2),
           l2S=mean(l2S), l2S_sd=sd(l2S),
           lF_Omega=mean(lF_Omega), lF_Omega_sd=sd(lF_Omega),
           fpr_opt=mean(fpr_opt), fpr_opt_sd=sd(fpr_opt),
           tpr_opt=mean(tpr_opt), tpr_opt_sd=sd(tpr_opt),
           fpr=colMeans(fpr), tpr=colMeans(tpr))
  res <- round(res, 3)
  
  return (res)
}

run_iter <- function(data, estimator, nlambda, REML, eta, eta_gs, maxits,
                     standardize, progress){
  X_hat <- get_X_hat(data$data, estimator=estimator)
  dat <- data.frame(y=data$data$y,
                    X=cbind(data$data$t, X_hat),
                    cluster=data$data$cluster)
  
  ic.fit <- ic.lmm_lasso(data$formula, dat,
                         p=data$data$p, q=data$data$q, d=data$data$d,
                         n=data$data$n, N=data$data$N,
                         theta_0=data$init[[estimator]]$theta_0,
                         beta_0=data$init[[estimator]]$beta_0,
                         sig2e_0=data$init[[estimator]]$sig2e_0,
                         lambda.max=data$lambda.max,
                         nlambda=nlambda,
                         penalty.factor=c(0, rep(1, data$data$p_x)),
                         REML=REML, eta=eta, eta_gs=eta_gs,
                         maxits=maxits, standardize=standardize,
                         progress=progress)
  
  fpr <- rowMeans(ic.fit$S_hat[, -data$data$S])
  tpr <- rowMeans(ic.fit$S_hat[, data$data$S])
  
  fpr_opt <- mean(ic.fit$beta_opt[-data$data$S] != 0)
  tpr_opt <- mean(ic.fit$beta_opt[data$data$S] != 0)
  
  l2 <- norm( ic.fit$beta_opt - data$data$beta, "2" )^2
  l2S <- norm( ic.fit$beta_opt[data$data$S] - data$data$beta[data$data$S], "2" )^2
  l2_mu <- norm( ic.fit$mu_opt - data$data$mu, "2" )^2
  
  Omega_0_hat <- ic.fit$Omega_0_opt - get_OmegaBias(data$data, ic.fit$beta_opt,
                                                     estimator)
  # get nearest positive definite matrix:
  Omega_0_hat <- as.matrix(Matrix::nearPD(Omega_0_hat)$mat)
  
  lF_Omega <- norm( data$data$Omega_0 - Omega_0_hat, "F" )^2
  
  res <- c(l2_mu=l2_mu, l2=l2, l2S=l2S, lF_Omega=lF_Omega,
           fpr_opt=fpr_opt, tpr_opt=tpr_opt,
           fpr=fpr, tpr=tpr)
  
  return ( res )
}