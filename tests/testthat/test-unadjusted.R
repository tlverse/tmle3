context("Basic interventions: TSM for single static intervention")

library(sl3)
# library(tmle3)
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
  W = c(),
  A = "parity01",
  Y = "haz01"
)

Q_learner <- make_learner(Lrnr_glm)
g_learner <- make_learner(Lrnr_mean)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")
targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

# define parameter
intervention <- define_lf(LF_static, "A", value = 1)
tsm <- define_param(Param_TSM, targeted_likelihood, intervention)
updater$tmle_params <- tsm
# define update method (submodel + loss function)


# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(tsm), updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# task for A=1
cf_task <- tsm$cf_likelihood$cf_tasks[[1]]

# get Q
EY1 <- likelihood$get_likelihoods(cf_task, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task, "Y")
EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

# get G
pA1 <- likelihood$get_likelihoods(cf_task, "A")
pDelta1 <- cbind(pA1, pA1)

W <- 0 * Q # need something here so tmle doesn't break but it shouldn't be used
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = NULL,
  W = W,
  Delta = tmle_task$get_tmle_node("A"),
  Q = Q,
  pDelta1 = pDelta1,
  family = "binomial",
  alpha = 0.995,
  target.gwt = FALSE,
  prescreenW.g = FALSE
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
