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
source("simfunctions.R")
save_dir <- getwd()
output_file <- "progress_main.txt"
cat("", file=output_file)

set.seed(281)

load("data_pregen.RData")

opt_grid <- read.csv("options.csv", header=TRUE)

res_len <- nrow(opt_grid)
res_col <- ncol(opt_grid)

# currently this code works only if estimator is the only variable changing
results <- opt_grid
for (r in 1:res_len){
  opt <- opt_grid[r,]
  
  cat("condition", r, "\nn", opt$n, "| k", opt$k, "| tNSR", opt$tnsr,
      "| structure", opt$structure,
      "| estimator", opt$estimator, "\n", file=output_file, append=TRUE)
  
  res <- run_sim(data=data, estimator=opt$estimator, runs=length(data),
                 of=output_file, nlambda=20, REML=TRUE, eta=1e-4, eta_gs=5e-3,
                 maxits=5e4, standardize=TRUE, progress=FALSE, parallel=TRUE,
                 mc.cores=60)
  
  results[r, (res_col+1):(res_col+length(res))] <- res
  
  write.table(results[r,], file=paste(save_dir, "results.csv", sep="/"),
              sep=",", quote=FALSE, append=TRUE,
              col.names=FALSE, row.names=FALSE)
  
  cat("results updated\n\n", file=output_file, append=TRUE)
}