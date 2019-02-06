context("Direct ATE for single node interventions.")

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

logit_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial, loss_loglik_binomial)
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# LF_fit$undebug("get_likelihood")
# estimate likelihood
factor_list <- list(
  define_lf(LF_emp, "W", NA),
  define_lf(LF_fit, "A", type = "density", learner = learner_list[["A"]]),
  define_lf(LF_fit, "Y", type = "mean", learner = learner_list[["Y"]])
)

likelihood_def <- Likelihood$new(factor_list)
initial_likelihood <- likelihood_def$train(tmle_task)

updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define parameter
intervention_treatment <- define_lf(LF_static, "A", value = 1)
intervention_control <- define_lf(LF_static, "A", value = 0)
ate <- define_param(Param_ATE, targeted_likelihood, intervention_treatment, intervention_control)
updater$tmle_params <- list(ate)


# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(ate), updater)


# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
tmle3_epsilon <- updater$epsilons[[1]]$Y

submodel_data <- updater$generate_submodel_data(initial_likelihood, tmle_task, -1)
#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# tasks for A=1 and A=0
cf_task1 <- ate$cf_likelihood_treatment$cf_tasks[[1]]
cf_task0 <- ate$cf_likelihood_control$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task1, "Y")
EY0_final <- targeted_likelihood$get_likelihoods(cf_task0, "Y")
EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = cbind(tmle_task$get_tmle_node("W"),tmle_task$get_tmle_node("W")),
  Q = Q,
  g1W = pA1
)


# extract estimates
classic_psi <- tmle_classic_fit$estimates$ATE$psi
classic_se <- sqrt(tmle_classic_fit$estimates$ATE$var.psi)
tol <- 1/sqrt(tmle_task$nrow)
test_that("psi matches result from classic package", expect_equal(tmle3_psi, classic_psi, tol=tol))
test_that("se matches result from classic package", expect_equal(tmle3_se, classic_se, tol=tol))
