context("Basic interventions: TSM for single static intervention")

library(sl3)
# library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)
# setup data for test
set.seed(1234)
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)

node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  A = "parity01",
  Y = "haz01"
)

# drop missing A for now, might add back to test later
missing_W <- apply(is.na(data[, c(node_list$W, node_list$A), with = FALSE]), 1, any)
data <- data[!missing_W]

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
learner_list <- list(Y = Q_learner, A = g_learner, delta_Y = Q_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
Q_task <- tmle_task$get_regression_task("Y", drop_censored = TRUE)
Q_learner <- learner_list$Y
# debug_train(Q_learner)
Q_fit <- Q_learner$train(Q_task)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = FALSE)

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
intervention1 <- define_lf(LF_static, "A", value = 1)
intervention0 <- define_lf(LF_static, "A", value = 0)
tsm1 <- define_param(Param_TSM, targeted_likelihood, intervention1)
tsm0 <- define_param(Param_TSM, targeted_likelihood, intervention0)
ate <- define_param(Param_delta, targeted_likelihood, delta_param_ATE, list(tsm0, tsm1))
params <- list(tsm0, tsm1, ate)
updater$tmle_params <- params
H0W <- tsm0$clever_covariates(tmle_task)$Y
H1W <- tsm1$clever_covariates(tmle_task)$Y
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, params, updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est[3]
tmle3_se <- tmle_fit$summary$se[3]
tmle3_epsilon <- updater$epsilons[[1]]$Y

submodel_data <- updater$generate_submodel_data(
  initial_likelihood, tmle_task,
  "full"
)

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

cf_task1 <- tsm1$cf_likelihood$cf_tasks[[1]]
cf_task0 <- tsm0$cf_likelihood$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task1, "Y")
EY0_final <- targeted_likelihood$get_likelihoods(cf_task0, "Y")
Qstar <- cbind(EY0_final, EY1_final)
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
Q <- cbind(EY0, EY1)
Delta <- tmle_task$get_tmle_node("delta_Y")

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")
pDelta1a1 <- initial_likelihood$get_likelihoods(cf_task1, "delta_Y")
pDelta1a0 <- initial_likelihood$get_likelihoods(cf_task0, "delta_Y")

g10W <- pA1 * (pDelta1a1)
g00W <- (1 - pA1) * (pDelta1a0)

pDelta1 <- cbind(pDelta1a0, pDelta1a1)
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = tmle_task$get_tmle_node("W"),
  Delta = Delta,
  Q = Q,
  g1W = pA1,
  pDelta1 = pDelta1,
  family = "binomial",
  alpha = 0.995,
  target.gwt = FALSE,
  gbound=0.025
)


# extract estimates
classic_psi <- tmle_classic_fit$estimates$ATE$psi
classic_se <- sqrt(tmle_classic_fit$estimates$ATE$var.psi)
classic_epsilon <- tmle_classic_fit$epsilon
classic_Qstar <- tmle_classic_fit$Qstar

test_that("Qstar matches result from classic package", {
  expect_equivalent(Qstar, classic_Qstar)
})

test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi)
})

test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se)
})

test_that("epsilon matches resullt from classic package", {
  expect_equivalent(tmle3_epsilon, classic_epsilon)
})
