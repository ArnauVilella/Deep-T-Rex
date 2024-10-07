#!/usr/bin/env Rscript

library(optparse)
library(progress)
library(TRexSelector)

option_list <- list(
  make_option(c("--target_FDR"), type="numeric", default=0.1, 
              help="Target FDR [default %default]"),
  make_option(c("--T_stop"), type="integer", default=1, 
              help="T_stop [default %default]"),
  make_option(c("--N_data"), type="integer", default=100, 
              help="Dataset size to create [default %default]"),
  make_option(c("--n"), type="integer", default=75, 
              help="Number of observations [default %default]"),
  make_option(c("--p"), type="integer", default=150, 
              help="Number of variables [default %default]"),
  make_option(c("--K"), type="integer", default=100, 
              help="Number of random experiments [default %default]"),
  make_option(c("--num_act"), type="integer", default=3, 
              help="Number of active variables [default %default]"),
  make_option(c("--num_dummies"), type="integer", default=100, 
              help="Number of dummies [default %default]"),
  make_option(c("--SNR"), type="numeric", default=1.0, 
              help="SNR [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

target_FDR <- opt$target_FDR
T_stop <- opt$T_stop
N_data <- opt$N_data
n <- opt$n
p <- opt$p
K <- opt$K
num_act <- opt$num_act
num_dummies <- opt$num_dummies
SNR <- opt$SNR

create_subdir <- function(subdir) {
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
  }
}

FDP_inv <- function(target_FDR, phi_T, true_actives, K, T_stop=T_stop) {
  min_v = 1
  for (v in seq(1, 0.5, by=-1/K)) {
    FDP <- sum(phi_T[-true_actives, T_stop] > v) / max(sum(phi_T[, T_stop] > v), 1)
    if (FDP < target_FDR){
      min_v = v
    }
  }
  return(min_v)
}

main_dir <- "data"
create_subdir(main_dir)

pb <- progress_bar$new(format = "Generating data [:bar] :current/:total (:percent)", total=N_data)

for (i in 1:N_data) {
  beta <- sample(c(rep(1, times=num_act), rep(0, times=p - num_act)))
  true_actives <- which(beta > 0)
  
  X <- matrix(stats::rnorm(n * p), nrow=n, ncol=p)
  sd <- sqrt(var(X %*% beta)/SNR)
  y <- X %*% beta + stats::rnorm(n, sd=sd)
  
  res_exp <- random_experiments(X, y, K=K)
  Phi_mat <- res_exp$phi_T_mat  # It's transposed, so T_stop columns and p rows
  v_thresh <- FDP_inv(target_FDR, res_exp$phi_T_mat, true_actives, K)
  
  subdirs <- list(X="X", y="y", beta="beta", true_actives="true_actives", 
                  Phi_mat="Phi_mat", v_thresh="v_thresh")
  
  for (name in names(subdirs)) {
    subdir <- file.path(main_dir, subdirs[[name]])
    create_subdir(subdir)
    
    file_name <- paste0(subdir, "/", name, "_", i, ".txt")
    
    write.table(get(name), file = file_name, row.names=FALSE, col.names=FALSE)
  }
  
  pb$tick(1)
}