context("Stratification - estimate TSM in strata")

library(sl3)
library(uuid)
library(assertthat)
library(data.table)
library(future)
# data
washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv", stringsAsFactors = TRUE)
node_list <- list(
  W = c(
    "month", "aged", "momage",
    "momheight", "hfiacat"
  ),
  V = "Nlt18",
  A = "tr",
  Y = "whz"
)
processed <- process_missing(washb_data[1:500,], node_list)
data <- processed$data
node_list <- processed$node_list
# leaners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_xgboost <- make_learner(Lrnr_xgboost)

ls_metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_xgboost),
  metalearner = ls_metalearner
)
sl_A <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_xgboost),
  metalearner = mn_metalearner
)

learner_list <- list(A = sl_A, Y = sl_Y)
# estimators
tmle_spec <- tmle_MSM()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new()

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_param <- tmle_spec$make_params(tmle_task, targeted_likelihood)
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_param, updater)

tmle_ests <- tmle_fit$summary$tmle_est
