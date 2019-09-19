context("Repeated measures: check variance is computed from subject-level ICs")

library(sl3)
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

Q_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

g_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, Q_lib, metalearner)
g_learner <- make_learner(Lrnr_sl, g_lib, metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list, id = "subjid")

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define parameter
# cf_likelihood <- CF_Likelihood$new(likelihood, intervention)

# define update method (submodel + loss function)
# updater <- tmle_spec$make_updater(likelihood, list(tsm))
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
intervention <- define_lf(LF_static, "A", value = 1)
intervention_2 <- define_lf(LF_static, "A", value = 0)

tsm <- define_param(Param_TSM, targeted_likelihood, intervention)
tsm_2 <- define_param(Param_TSM, targeted_likelihood, intervention_2)

updater$tmle_params <- list(tsm)

tmle_fit <- fit_tmle3(
  tmle_task, targeted_likelihood, updater$tmle_params,
  updater
)

# extract results
# debugonce(summary_from_estimates)
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# task for A=1
# cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(),
# data.table(A = 1))
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
  pDelta1 = pDelta1,
  id = tmle_task$id
)

# extract estimates
classic_psi <- tmle_classic_fit$estimates$EY1$psi
classic_se <- sqrt(tmle_classic_fit$estimates$EY1$var.psi)

# only approximately equal (although it's O(1/n))
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol = 1e-3)
})

# only approximately equal (although it's O(1/n))
test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol = 1e-3)
})
