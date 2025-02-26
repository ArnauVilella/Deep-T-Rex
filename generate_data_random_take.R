#!/usr/bin/env Rscript

library(optparse)
library(progress)
library(TRexSelector)
options(scipen=999)

option_list <- list(
  make_option(c("--target_FDR"), type="numeric", default=0.1, 
              help="Target FDR [default %default]"),
  make_option(c("--T_stop"), type="integer", default=1, 
              help="T_stop [default %default]"),
  make_option(c("--N_systems"), type="integer", default=100, 
              help="Number of different dystems that create the dataset [default %default]"),
  make_option(c("--n"), type="integer", default=75, 
              help="Number of observations [default %default]"),
  make_option(c("--p"), type="integer", default=150, 
              help="Number of variables [default %default]"),
  make_option(c("--K"), type="integer", default=100, 
              help="Number of random experiments [default %default]"),
  make_option(c("--num_act"), type="integer", default=3, 
              help="Number of active variables [default %default]"),
  make_option(c("--num_dummies"), type="integer", default=150, 
              help="Number of dummies [default %default]"),
  make_option(c("--SNR"), type="numeric", default=1.0, 
              help="SNR [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

target_FDR <- opt$target_FDR
T_stop <- opt$T_stop
N_systems <- opt$N_systems
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

main_dir <- "data"
create_subdir(main_dir)

pb <- progress_bar$new(format = "Generating data [:bar] :current/:total (:percent)", total=N_systems)

for (i in 1:N_systems) {
  beta <- sample(c(rep(1, times=num_act), rep(0, times=p - num_act)))
  true_actives <- which(beta > 0)
  
  X <- matrix(stats::rnorm(n * p), nrow=n, ncol=p)
  sd <- sqrt(var(X %*% beta)/SNR)
  y <- X %*% beta + stats::rnorm(n, sd=sd)
  
  res_exp <- random_experiments(X, y, K=K, T_stop=T_stop)
  Phi_mat <- res_exp$phi_T_mat  # It's transposed, so T_stop columns and p rows
  Phi_mat <- Phi_mat[, ncol(Phi_mat)]
  v <- runif(1, min=0.5, max=1)
  FDR <- (sum((1 - beta) * Phi_mat > v) / max(1, sum(Phi_mat > v)))
  
  subdirs <- list(beta="beta", Phi_mat="Phi_mat", v="v", FDR="FDR")
  
  for (name in names(subdirs)) {
    subdir <- file.path(main_dir, subdirs[[name]])
    create_subdir(subdir)
    
    file_name <- paste0(subdir, "/", name, "_", i, ".txt")
    
    write.table(get(name), file = file_name, row.names=FALSE, col.names=FALSE)
  }
  
  pb$tick(1)
}