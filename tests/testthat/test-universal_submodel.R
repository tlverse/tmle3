context("Universal submodel")

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
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
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
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle3_Update$new(
  one_dimensional = TRUE, constrain_step = TRUE,
  maxit = 10000, cvtmle = TRUE,
  convergence_type = "sample_size"
)

# updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
intervention <- define_lf(LF_static, "A", value = 1)
# params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
tsm <- define_param(Param_TSM, targeted_likelihood, intervention)

updater$tmle_params <- tsm

tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tsm, updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
tmle3_epsilon <- updater$epsilons[[1]]$Y


#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# task for A=1
# cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A = 1))
cf_task <- tsm$cf_likelihood$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task, "Y")
EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task, "A")
pDelta1 <- cbind(pA1, pA1)
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = NULL,
  W = tmle_task$get_tmle_node("W"),
  Delta = tmle_task$get_tmle_node("A"),
  Q = Q,
  pDelta1 = pDelta1
)

cf_task <- tsm$cf_likelihood$cf_tasks[[1]]
# debugonce(tsm$cf_likelihood$enumerate_cf_tasks)

# extract estimates
classic_psi <- tmle_classic_fit$estimates$EY1$psi
classic_se <- sqrt(tmle_classic_fit$estimates$EY1$var.psi)
classic_epsilon <- tmle_classic_fit$epsilon[["H1W"]]
classic_Qstar <- tmle_classic_fit$Qstar[, 2]

# test_that("Qstar matches result from classic package", expect_equivalent(EY1_final, classic_Qstar))
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol = 2e-3)
})
test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol = 1e-3)
})
