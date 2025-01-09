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
g_treatment <- rep(0.5, n_obs)
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

# define an adaptive design for evaluation
g_adapt_1 <- g_treatment / 2 * (W > 0.5) + (1 - g_treatment / 2) * (W < 0.5)

# Test 1
## Define tmle_spec
tmle_spec <- tmle3_Spec_ADSM$new(
  treatment_level = 1,
  control_level = 0,
  g_treat = g_treatment,
  g_adapt = g_adapt_1
)

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
truth_1 <- mean(g_adapt_1 * EY1 + (1 - g_adapt_1) * EY0)

test_that("TMLE CI includes truth", {
  expect_lte(abs(truth_1 - tmle_fit$summary$tmle_est), tmle_fit$summary$se * 1.96)
})


# Test 2
## More extreme g_adapt
g_adapt_2 <- rep(0.99, n_obs)

tmle_spec_2 <- tmle_ADSM(
  treatment_level = 1,
  control_level = 0,
  g_treat = g_treatment,
  g_adapt = g_adapt_2
)

tmle_task_2 <- tmle_spec_2$make_tmle_task(data, node_list)
initial_likelihood_2 <- tmle_spec_2$make_initial_likelihood(
  tmle_task_2,
  learner_list
)
targeted_likelihood_2 <- Targeted_Likelihood$new(initial_likelihood_2)
tmle_params_2 <- tmle_spec_2$make_params(tmle_task_2, targeted_likelihood_2)
tmle_fit_2 <- fit_tmle3(
  tmle_task_2, targeted_likelihood_2, tmle_params_2,
  targeted_likelihood_2$updater
)

truth_2 <- mean(g_adapt_2 * EY1 + (1 - g_adapt_2) * EY0)

test_that("TMLE CI includes truth", {
  expect_lte(abs(truth_2 - tmle_fit_2$summary$tmle_est), tmle_fit_2$summary$se * 1.96)
})


# Test 3
## Compare ADSM estimate where p(A = 1) = 0.99 with TSM estimate for E[Y_{A=1}]
tsm_spec <- tmle_TSM_all()
tmle_task_tsm <- tmle_spec$make_tmle_task(data, node_list)
initial_likelihood_tsm <- tmle_spec$make_initial_likelihood(
  tmle_task_tsm,
  learner_list
)

targeted_likelihood_tsm <- Targeted_Likelihood$new(initial_likelihood_tsm)
tmle_params_tsm <- tsm_spec$make_params(tmle_task_tsm, targeted_likelihood_tsm)
tmle_fit_tsm <- fit_tmle3(
  tmle_task_tsm, targeted_likelihood_tsm, tmle_params_tsm,
  targeted_likelihood_tsm$updater
)
whTrt1 <- which(tmle_fit_tsm$tmle_param_names == "E[Y_{A=1}]")

test_that("TMLE Estimate almost equals TSM(A = 1)", {
  expect_lte(abs(tmle_fit_tsm$summary$tmle_est[whTrt1] - tmle_fit_2$summary$tmle_est), 0.01)
})

# Test 4
## Compare ADSM under ODTR: I(CATE > 0)
g_adapt_4 <- as.numeric(I(EY1 > EY0))

tmle_spec_4 <- tmle_ADSM(
  treatment_level = 1,
  control_level = 0,
  g_treat = g_treatment,
  g_adapt = g_adapt_4
)

tmle_task_4 <- tmle_spec_4$make_tmle_task(data, node_list)
initial_likelihood_4 <- tmle_spec_4$make_initial_likelihood(
  tmle_task_4,
  learner_list
)
targeted_likelihood_4 <- Targeted_Likelihood$new(initial_likelihood_4)
tmle_params_4 <- tmle_spec_4$make_params(tmle_task_4, targeted_likelihood_4)
tmle_fit_4 <- fit_tmle3(
  tmle_task_4, targeted_likelihood_4, tmle_params_4,
  targeted_likelihood_4$updater
)

truth_4 <- mean(g_adapt_4 * EY1 + (1 - g_adapt_4) * EY0)

test_that("TMLE CI includes truth", {
  expect_lte(abs(truth_4 - tmle_fit_4$summary$tmle_est), tmle_fit_4$summary$se * 1.96)
})

