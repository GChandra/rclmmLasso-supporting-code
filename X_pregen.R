library(Matrix)
library(MASS)
library(matrixcalc)
library(nloptr)
library(lme4)
library(tidyr)
library(dplyr)
library(parallel)


# Set to appropriate working directory
setwd("supporting_code")

source("cdlmmLasso.R")
source("helperfunctions.R")
source("meMethods.R")
save_dir <- getwd()

run_gen <- function(r, n, m, p1, p2, q, k, mu, beta_t, beta1, beta2,
                    mu_x, Sigma_x, Omega_0, Sigma_d1, z, Sigma_e, theta,
                    estimators){
  data <- data_gen_lmm_bsl(n=n, m=m, p1=p1, p2=p2, q=q, k=k, mu=mu,
                           beta_t=beta_t, beta1=beta1, beta2=beta2,
                           mu_x=mu_x, Sigma_x=Sigma_x, Omega_0=Omega_0,
                           Sigma_d1=Sigma_d1, z=z, Sigma_e=Sigma_e,
                           theta=theta)
  
  init <- list()
  lambda.ub <- -Inf
  for (estim in estimators){
    X_hat <- get_X_hat(data, estimator=estim)
    dat <- data.frame(y=data$y,
                      X=cbind(data$t, X_hat),
                      cluster=data$cluster)
    frm <- paste("y ~ 1 + ", paste( names(dat)[grepl("X.", names(dat))],
                                    collapse=" + " ),
                 " + (1 | cluster)",
                 sep="" )
    init[[estim]] <- ic.lmm_lasso(frm, dat,
                                  data$p, data$q, d=data$d, n=data$n, N=data$N,
                                  penalty.factor=c(0, rep(1, data$p_x)),
                                  REML=TRUE, eta=1e-4, eta_gs=5e-3, maxits=5e4,
                                  standardize=TRUE, getInit=TRUE)
    
    if (init[[estim]]$lambda.max > lambda.ub){
      lambda.ub <- init[[estim]]$lambda.max
    }
  }
  
  return(list(
    data=data,
    formula=frm,
    init=init,
    lambda.ub=lambda.ub
  ))
}

set.seed(281)

# values below chosen to reflect IDATA data (n may be different, though)
opt_grid <- read.csv("options.csv", header=TRUE)
estimators <- opt_grid$estimator

# I am assuming n, k, tnsr, structure don't change in the conditions.
opt <- opt_grid[1,]

n <- opt$n
m <- 3
k <- opt$k
p1 <- 30
p2 <- 2
q <- 1 # random effect dimension
s1 <- 5

beta_t <- 0

S1 <- sample(1:p1, size=s1)
beta1 <- rep(0, p1)
beta1[S1] <- rnorm(s1, mean=0.8, sd=1)
beta2 <- rnorm(2, mean=0.8, sd=1)

sig2e <- 0.6
Sigma_e <- sig2e * diag(m)

sig2b <- 15
Omega_0 <- matrix(c(sig2b),
                  nrow=q, ncol=q)
Omega_0_c <- t(chol(Omega_0)) / sqrt(sig2e)
theta <- Omega_0_c[lower.tri(Omega_0_c, diag=TRUE)]

z <- c(-1)

mu <- 20
mu_x <- c(rep(4, p1+p2))
mu_x[(p1+1):(p1+p2)] <- c(63, 0.5)

runs <- 500

sig2x1 <- 1
if (opt$structure == "diagonal"){
  Sigma_x1x1 <- choose_matrix(sig2=sig2x1, p=p1, structure="diag")
  Sigma_x1x2 <- matrix(0, nrow=p1, ncol=p2)
  Sigma_x2x2 <- matrix(c(35, 0, 0, 0.25), nrow=p2, ncol=p2)
} else if (opt$structure == "ar1"){
  Sigma_x1x1 <- choose_matrix(sig2=sig2x1, rho=0.65, p=p1, structure="AR1")
  Sigma_x1x2 <- matrix(0, nrow=p1, ncol=p2)
  Sigma_x1x2[1,2] <- -0.2*sqrt(0.25)
  Sigma_x2x2 <- matrix(c(35, 0, 0, 0.25), nrow=p2, ncol=p2)
} else {
  stop("unsupported matrix structure '", opt$structure, "'")
}
Sigma_x <- rbind(
  cbind(Sigma_x1x1, Sigma_x1x2),
  cbind(t(Sigma_x1x2), Sigma_x2x2)
)

sig2d1 <- opt$k * opt$tnsr / (1-opt$tnsr) * sig2x1
Sigma_d1 <- sig2d1 * diag(p1)

data <- mclapply(1:runs, run_gen, n=n, m=m, p1=p1, p2=p2, q=q, k=k,
                 mu=mu, beta_t=beta_t, beta1=beta1, beta2=beta2,
                 mu_x=mu_x, Sigma_x=Sigma_x, Omega_0=Omega_0, Sigma_d1=Sigma_d1,
                 z=z, Sigma_e=Sigma_e, theta=theta, estimators=estimators,
                 mc.cores=60)

lambda.max <- quantile( sapply(data, function(x) x$lambda.ub), probs=0.05 )
for (i in 1:length(data)){
  data[[i]][["lambda.max"]] <- lambda.max
}

save(data, file="data_pregen.RData")
