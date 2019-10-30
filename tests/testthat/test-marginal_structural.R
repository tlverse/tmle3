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
  V = "sex", #Nlt18
  A = "tr",
  Y = "whz"
)
processed <- process_missing(washb_data[1:500,], node_list)
data <- processed$data
node_list <- processed$node_list
# leaners
qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_xgboost"
)
glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_xgboost"
)
ls_metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, ls_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
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
pA <- 1/tmle_fit$tmle_params[[2]]$strata$weight
wm <- weighted.mean(tmle_ests[-1],pA)
test_that("overall ATE is weighted average of strata ATEs",expect_equal(tmle_ests[[1]],wm))

ses <- tmle_fit$summary$se

test_that("overall ATE has lower SE than strata ATEs",expect_equal(which.min(ses),1))
