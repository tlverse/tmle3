context("Odds Ratios")

library(sl3)
library(uuid)
library(assertthat)
library(data.table)
library(future)
set.seed(1234)

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

# metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_glm_fast)
g_learner <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_OR(baseline = 0, contrast = 1)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define param
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)

# extract results
summary <- tmle_fit$summary
tmle3_psi <- tmle_fit$summary$psi_transformed[[3]]
tmle3_se <- tmle_fit$summary$se[[3]]
#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# tasks for A=1 and A=0
cf_task1 <- tmle_params[[2]]$cf_likelihood$cf_tasks[[1]]
cf_task0 <- tmle_params[[1]]$cf_likelihood$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task1, "Y")
EY0_final <- targeted_likelihood$get_likelihoods(cf_task0, "Y")
Q <- cbind(EY0, EY1)

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = cbind(tmle_task$get_tmle_node("W"), tmle_task$get_tmle_node("W")),
  Q = Q,
  g1W = pA1,
  family = "binomial"
)


# extract estimates
classic_psi <- tmle_classic_fit$estimates$OR$psi
classic_se <- sqrt(tmle_classic_fit$estimates$OR$var.log.psi)
tol <- 1 / sqrt(tmle_task$nrow)
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol = tol)
})

test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol = tol)
})
