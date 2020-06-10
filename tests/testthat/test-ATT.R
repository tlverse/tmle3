context("ATT interventions: treatment effect amongst the treated")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# setup data for test
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0
node_list <- list(
  W = c("sexn"),
  A = "parity01",
  Y = "haz01"
)

qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_ATT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# LF_fit$undebug("get_likelihood")
# estimate likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

updater <- tmle3_Update$new(
  cvtmle = FALSE, convergence_type = "sample_size",
  constrain_step = TRUE, one_dimensional = TRUE, delta_epsilon = 0.001,
  optim_delta_epsilon = FALSE
)

# debugonce(updater$update_step)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)


# define parameter
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
att <- tmle_params[[1]]

# fit tmle update
tmle_fit <- fit_tmle3(
  tmle_task, targeted_likelihood, list(att), updater,
  max_it
)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
tmle3_epsilon <- updater$epsilons[[1]]$Y

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# tasks for A=1 and A=0
cf_task1 <- att$cf_likelihood_treatment$cf_tasks[[1]]
cf_task0 <- att$cf_likelihood_control$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task1, "Y")
EY0_final <- targeted_likelihood$get_likelihoods(cf_task0, "Y")
# EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")

# debugonce(oneStepATT)
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = cbind(tmle_task$get_tmle_node("W"), tmle_task$get_tmle_node("W")),
  Q = Q,
  g1W = pA1,
  family = "binomial",
  alpha = 0.995,
  target.gwt = FALSE
)


# extract estimates
classic_psi <- tmle_classic_fit$estimates$ATT$psi
classic_se <- sqrt(tmle_classic_fit$estimates$ATT$var.psi)
tol <- 1 / sqrt(tmle_task$nrow)
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol)
})

test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol)
})
