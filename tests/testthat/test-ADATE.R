context("Test the TML Estimator for the Mean Outcome under a Counterfactual Adaptive Design")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)

set.seed(1234)

## simulate simple data for TML for adaptive design
n_obs <- 5000 # number of observations

## baseline covariates -- simple, binary
W <- rnorm(n_obs, 0, 1)

## create treatment based on baseline W
g_fair <- rep(0.5, n_obs)
g_treatment <- g_fair / 2 * (W > 0.5) + (1 - g_fair / 2) * (W < 0.5)
A <- sapply(g_treatment, rbinom, n = 1, size = 1)

EY1 <- W + W^2
EY0 <- W
Y <- A*EY1 + (1-A)*EY0 + rnorm(n_obs, mean = 0, sd = 1)

## organize data and nodes for tmle3
data <- data.table(W, A, Y)
node_list <- list(W = "W", A = "A", Y = "Y")

# learners used for conditional expectation regression (e.g., outcome)
mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
sl_lrnr <- Lrnr_sl$new(
  learners = list(mean_lrnr, glm_lrnr),
  metalearner = Lrnr_nnls$new()
)
learner_list <- list(A = mean_lrnr, Y = sl_lrnr)

# Test 1
## Define tmle_spec
tmle_spec <- tmle3_Spec_ADATE$new(
  treatment_level = 1,
  control_level = 0,
  g_treat = g_treatment)

## Define tmle task
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

## Make initial likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

## Create targeted_likelihood object
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)

## Define tmle param
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)

## Run TMLE
tmle_fit <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)
tmle_fit

## Truth
truth_1 <- mean(EY1-EY0)

test_that("TMLE CI includes truth", {
  expect_lte(abs(truth_1 - tmle_fit$summary$tmle_est), tmle_fit$summary$se * 1.96)
})

