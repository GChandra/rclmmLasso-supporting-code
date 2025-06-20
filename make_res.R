library(ggplot2)
library(tidyverse)

setwd("/Users/gchandra/Desktop/Ganesh/Research/project2/simulations/2025_06_04")

suffix <- ""

res <- read.csv(paste0("results", suffix, ".csv"), header=FALSE)

res_table <- res %>% reframe(
  n=V1,
  k=V2,
  tNSR=V3,
  structure=V4,
  estimator=V5,
  l2_mu = paste(V6, sprintf("(%.2f)", V7)),
  l2_beta = paste(V8, sprintf("(%.2f)", V9)),
  l2_beta_S = paste(V10, sprintf("(%.2f)", V11)),
  lF_Omega = paste(V12, sprintf("(%.2f)", V13)),
  fpr = paste(V14, sprintf("(%.2f)", V15)),
  tpr = paste(V16, sprintf("(%.2f)", V17))
)

write.csv(res_table, paste0("res_table", suffix, ".csv"))

nlambda <- (length(res) - 17) / 2
roc_res <- res %>%
  dplyr::select(estimator=V5, fpr.=18:(18+nlambda-1), tpr.=(18+nlambda):(18+2*nlambda-1))
roc_res <- roc_res %>%
  reshape(direction='long', 
          varying=matrix(
            colnames(dplyr::select(roc_res, contains("."))),
            nrow=2,
            byrow=TRUE
            ), 
          timevar='lambda',
          times=1:nlambda,
          v.names=c('fpr', 'tpr'))

ggplot(roc_res, aes(x=fpr, y=tpr, linetype=estimator, shape=estimator)) +
  geom_line() +
  geom_point(size=2) +
  labs(title="n=1000, p=33", x="False Positive Rate", y="True Positive Rate") +
  scale_linetype_manual(labels=c("_naive" = "Naive", "_rck"="RCk", "_rcu"="RCu",
                                 "_true"="True"),
                        values=c("solid", "dashed", "dotted", "dotdash"),
                        name="Estimator") +
  scale_shape_manual(labels=c("_naive" = "Naive", "_rck"="RCk", "_rcu"="RCu",
                                 "_true"="True"),
                        values=c(16, 15, 17, 3),
                        name="Estimator") +
  theme_bw()
