context("Different covariates: Specify adjustment sets in learners")

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

g_covariates <- "mage"
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner, covariates = g_covariates)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# verify that learners have expected covariate sets
g_fit <- initial_likelihood$factor_list[["A"]]$learner

test_that("g fit has a reduced covariate set", {
  expect_equal(g_fit$training_task$nodes$covariates, g_covariates)
  internal_glm_fit <- g_fit$fit_object$full_fit$fit_object$learner_fits[[1]]$fit_object$learner_fits[[2]]
  expect_equal(names(coef(internal_glm_fit)), c("intercept", g_covariates))
})

test_that("a manually constructed g fit matches", {
  g_task_manual <- sl3_Task$new(data, outcome = node_list$A, covariates = g_covariates, folds = tmle_task$folds)
  g_fit_manual <- learner_list$Y$train(g_task_manual)
  manual_preds <- g_fit_manual$predict()
  tmle_preds <- g_fit$predict()
  expect_equal(manual_preds, tmle_preds)
})

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = FALSE)

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
intervention <- define_lf(LF_static, "A", value = 1)
tsm <- define_param(Param_TSM, targeted_likelihood, intervention)
updater$tmle_params <- tsm

# debugonce(targeted_likelihood$update)
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(tsm), updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
tmle3_epsilon <- updater$epsilons[[1]]$Y

submodel_data <- updater$generate_submodel_data(
  initial_likelihood, tmle_task,
  "full"
)
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
  pDelta1 = pDelta1,
  family = "binomial",
  alpha = 0.995,
  target.gwt = FALSE
)

cf_task <- tsm$cf_likelihood$cf_tasks[[1]]
# debugonce(tsm$cf_likelihood$enumerate_cf_tasks)

# extract estimates
classic_psi <- tmle_classic_fit$estimates$EY1$psi
classic_se <- sqrt(tmle_classic_fit$estimates$EY1$var.psi)
classic_epsilon <- tmle_classic_fit$epsilon[["H1W"]]
classic_Qstar <- tmle_classic_fit$Qstar[, 2]

test_that("Qstar matches result from classic package", {
  expect_equivalent(EY1_final, classic_Qstar)
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
