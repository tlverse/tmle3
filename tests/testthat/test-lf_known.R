context("test-lf_known - Known likelihoods")
library(data.table)
library(sl3)
library(tmle3shift)
library(tmle3)
set.seed(12345)
sim_data <- function(n_obs = 1e3, n_w = 1, tx_mult = 2) {
  # n_obs - number of observations
  # n_w - number of baseline covariates
  # tx_mult - multiplier for the effect of W = 1 on the treatment
  
  # baseline covariates -- simple, binary
  W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))
  
  # create treatment based on baseline W
  A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))
  
  # create outcome as a linear function of A, W + white noise
  Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)
  
  # make output data
  data_obs <- data.table(W, A, Y)
  node_list <- list(W = "W", A = "A", Y = "Y")
  out <- list(data = data_obs, nodes = node_list)
  return(out)
}


sim_obj <- sim_data(1e6)
node_list <- sim_obj$nodes

tmle_spec <- tmle_shift(shift_val = 0.5)
tmle_task <- tmle_spec$make_tmle_task(sim_obj$data, sim_obj$nodes)

Q_task <- tmle_task$get_regression_task("Y")
g_task <- tmle_task$get_regression_task("A")

g_mean <- function(task){
  W <- g_task$data$W
  return(2*W)
}

g_dens <- function(task){
  mean_val <- g_mean(task)
  dens <- dnorm(task$Y, mean = mean_val)
}

Q_mean <- function(task){
  W <- task$data$W
  A <- task$data$A
  return(A+W)
}

# use known likelihoods
factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_known, "A", mean_fun = g_mean, density_fun = g_dens),
  define_lf(LF_known, "Y", mean_fun = Q_mean, type = "mean")
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)
tmle_params <- tmle_spec$make_params(tmle_task, likelihood)
estimates <- tmle_params[[1]]$estimates(tmle_task)
psi <- estimates$psi
var_EIF <- var(estimates$IC)
