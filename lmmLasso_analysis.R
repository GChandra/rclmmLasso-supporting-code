library(Matrix)
library(MASS)
library(matrixcalc)
library(nloptr)
library(lme4)
library(tidyverse)
library(pheatmap)
library(ggplot2)

select <- dplyr::select

setwd("supporting_code")
source("cdlmmLasso.R")

################################################################################
# load data
################################################################################
data <- read.csv("log_analysis_data.csv")
nutr_names <- names(data)[grepl("tn_", names(data))]
nutr_exc_names <- numeric(0)
nutr_inc_names <- setdiff(nutr_names, nutr_exc_names)

data <- data %>% select(-all_of(nutr_exc_names))

################################################################################
# estimate covariance matrices
################################################################################
# get nutrients at baseline (by ignoring duplicates)
bsl_nutr <- data[!duplicated(select(data, -c(bmi, month))), ] %>%
  select(-c(sex, task, bmi, month, age))

# split data frame by iid into list
bsl_list <- list()
for (id in unique(bsl_nutr$iid)){
  bsl_list[[id]] <- as.matrix(
    bsl_nutr[bsl_nutr$iid==id,] %>%
      select(-iid)
  )
}

n <- length(bsl_list)
p1 <- ncol(bsl_nutr) - 1
p2 <- 2 # age and se
k <- sapply(bsl_list, nrow, simplify=TRUE)

q <- 1 # number of random effects
d <- q*(q+1)/2 # number of variance components

Sigma_u_hat <- lapply(
  bsl_list,
  FUN = function(X) t(scale(X, scale=FALSE)) %*% scale(X, scale=FALSE)
)

Sigma_u_hat <- Reduce("+", Sigma_u_hat) / sum(k-1)

mu_x_hat <- mu_w_hat <- colMeans( select(bsl_nutr, -iid) )

Z <- data %>%
  select(iid, age, sex)
Z <- Z[!duplicated(Z), ]
mu_z_hat <- colMeans( select(Z, -iid) )

nu <- sum(k) - (sum(k^2) / sum(k))

Sigma_z_hat <- cov( select(Z, -iid) )

Sigma_xz_hat <- matrix(0, nrow=ncol(bsl_nutr)-1, ncol=p2)
Sigma_x_hat <- matrix(0, nrow=p1, ncol=p1)
Z_list <- list()
for (id in unique(Z$iid)){
  Z_list[[id]] <- as.matrix(
    Z[Z$iid==id,] %>%
      select(-iid)
  )
  
  Sigma_xz_hat <- Sigma_xz_hat +
    k[id] * matrix( colMeans( bsl_list[[id]] ) - mu_w_hat, ncol=1 ) %*%
    matrix( Z_list[[id]] - mu_z_hat, nrow=1 )
  
  Sigma_x_hat <- Sigma_x_hat +
    k[id] * matrix( colMeans( bsl_list[[id]] ) - mu_w_hat, ncol=1 ) %*%
    matrix( colMeans( bsl_list[[id]] ) - mu_w_hat, nrow=1 )
}
Sigma_xz_hat <- Sigma_xz_hat / nu
Sigma_x_hat <- ( Sigma_x_hat - (n-1)*Sigma_u_hat ) / nu

pheatmap(cov2cor(Sigma_x_hat), cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames=TRUE, show_colnames=TRUE)

pheatmap(cov2cor(Sigma_u_hat), cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames=TRUE, show_colnames=TRUE)

# noise-to-signal ratio based on diagonals
tnsr <- diag(Sigma_u_hat) / diag(Sigma_x_hat + Sigma_u_hat)

################################################################################
# LMM Lasso naive analysis
################################################################################
mn_data <- data %>%
  group_by(iid, bmi, age, sex, month) %>%
  summarise(across( .cols=ends_with("_asa24"),
                    .fns=mean,
                    .names="{.col}")
  )

frm <- paste("bmi ~ 1 + month + ", paste( names(mn_data)[grepl("tn_", names(mn_data))],
                                  collapse=" + "),
             " + age + sex + (1 | iid)",
             sep="" )

formula <- formula(frm)

frm <- lFormula(eval(formula), mn_data)
Z <- t(frm$reTrms$Zt)

intercept <- all(frm$X[,1]==1)
y <- frm$fr[,1]
X <- frm$X[,-1]

penalty.factor <- rep(1, ncol(X))
penalty.factor[1] <- 0

fit.ic.naive <- ic.lmm_lasso(
  formula, mn_data,
  p=ncol(X), q=q, d=d, n=length(unique(mn_data$iid)), N=nrow(X),
  penalty.factor=penalty.factor,
  REML=TRUE, BIC=TRUE,
  standardize=TRUE, progress=2, getInit=FALSE,
  control=lmerControl(optCtrl=list(xtol_rel=1e-5))
)

# univariate lmer analysis Naive for all nutrients
fit.naive.uni <- list()
for (nutr in nutr_inc_names){
  fit.naive.uni[[nutr]] <- lmer(paste("bmi ~", nutr, "+ (1 | iid)"), data=mn_data)
}

# lmer analysis for Naive (i.e. with no penalty)
fit.lmer.naive <- lmer(
  frm$formula, mn_data,
  verbose=0,
  control=lmerControl(optCtrl=list(xtol_rel=1e-5))
)
################################################################################
# LMM Lasso RCu analysis
################################################################################
min_task <- data %>%
  group_by(iid) %>%
  summarise(min_task=min(task))

rcu_data <- merge(data, min_task, by="iid") %>%
  filter(task==min_task) %>% select(-c(task, min_task))

A <- rbind(Sigma_x_hat,
           t(Sigma_xz_hat)
)

B <- rbind(
  cbind(Sigma_x_hat + Sigma_u_hat, Sigma_xz_hat),
  cbind(t(Sigma_xz_hat), Sigma_z_hat)
)
Lambda_hat <- solve(B, A)

pheatmap(abs(Lambda_hat), cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames=TRUE, show_colnames=TRUE)

pheatmap(Lambda_hat, cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames=TRUE, show_colnames=TRUE)

rcu_data[, row.names(Sigma_x_hat)] <-
  rep(1,nrow(rcu_data)) %*% (t(mu_w_hat) - ( cbind(t(mu_w_hat), t(mu_z_hat)) %*% Lambda_hat )) +
  as.matrix(rcu_data[, c(row.names(Sigma_x_hat), "age", "sex")]) %*% Lambda_hat

frm <- frm <- paste("bmi ~ 1 + month + ", paste( names(rcu_data)[grepl("tn_", names(rcu_data))],
                                         collapse=" + "),
                    " + age + sex + (1 | iid)",
                    sep="" )

formula <- formula(frm)

frm <- lFormula(eval(formula), rcu_data)
Z <- t(frm$reTrms$Zt)

intercept <- all(frm$X[,1]==1)
y <- frm$fr[,1]
X <- frm$X[,-1]

penalty.factor <- rep(1, ncol(X))
penalty.factor[1] <- 0

fit.ic.rcu <- ic.lmm_lasso(
  formula, rcu_data,
  p=ncol(X), q=q, d=d, n=length(unique(rcu_data$iid)), N=nrow(X),
  lambda.min.ratio=0.01,
  penalty.factor=penalty.factor,
  REML=TRUE, BIC=TRUE,
  standardize=TRUE, progress=2, getInit=FALSE,
  control=lmerControl(optCtrl=list(xtol_rel=1e-5))
)

fit.lmer.rcu <- lmer(
  frm$formula, rcu_data,
  verbose=1,
  control=lmerControl(optCtrl=list(xtol_rel=1e-5))
)
################################################################################
# put results into tables
################################################################################
beta_hat.lasso.naive <- fit.ic.naive$beta_opt
beta_hat.lasso.rcu <-  fit.ic.rcu$beta_opt
lmmLasso.beta_hat <- data.frame(naive=round(beta_hat.lasso.naive,3),
                                RCu=round(beta_hat.lasso.rcu,3)) %>%
  arrange(desc(abs(naive)), desc(abs(RCu)))

beta_hat.lmm.naive <- fixef(fit.lmer.naive)[-1]
beta_hat.lmm.rcu <- fixef(fit.lmer.rcu)[-1]
lmm.beta_hat <- data.frame(naive=beta_hat.lmm.naive,
                           RCu=beta_hat.lmm.rcu) %>%
  arrange(desc(abs(naive)), desc(abs(RCu)))

sel_nutrs <- which((fit.ic.rcu$beta_opt!=0) | (fit.ic.naive$beta_opt!=0))
sel_nutrs <- sel_nutrs[(sel_nutrs) > 1 & (sel_nutrs <= p1)] - 1
pheatmap(Lambda_hat[sel_nutrs,],
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames=TRUE, show_colnames=TRUE)

Lambda_hat_f <- cbind(
  Lambda_hat,
  solve(B, rbind(Sigma_xz_hat, Sigma_z_hat))
)
sig2_tau <- t(beta_hat.lasso.rcu[-1]) %*%
  B %*% ( diag(ncol(Lambda_hat_f)) - Lambda_hat_f ) %*%
  beta_hat.lasso.rcu[-1]
Omega_0_RC_hat <- fit.ic.rcu$Omega_0_opt - sig2_tau