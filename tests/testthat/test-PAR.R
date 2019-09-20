context("Population attributable risk and fraction")

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

q_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

g_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, q_lib, metalearner)
g_learner <- make_learner(Lrnr_sl, g_lib, metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_PAR(baseline_level = 1)

set.seed(1234)
# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle_spec$make_updater(convergence_type = "sample_size")

# define targeted_likelihood
targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

# define param
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
tmle_params[[1]]$estimates(tmle_task)
tmle_params[[1]]
# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)

# extract results
summary <- tmle_fit$summary

set.seed(1234)
data2 <- data.table::copy(data) # for data.table weirdness
spec <- tmle_PAR(baseline_level = 1)

tmle_fit_from_spec <- tmle3(
  spec, data2, node_list,
  learner_list
)

spec_summary <- tmle_fit_from_spec$summary

test_that("PAR manually and from Spec return the same results", {
  expect_equal(summary, spec_summary, tol = 1e-3)
})
