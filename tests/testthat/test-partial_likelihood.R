context("Fit subset of specified likelihood factors")

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

# define learners
Q_learner <- make_learner(Lrnr_glm_fast)
g_learner <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y = Q_learner, A = g_learner)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# test: fit likelihood values should only be available for factor subsets
test_that("Fitted likelihood values only available for W when subset is W", {
  tmle_spec <- tmle_RR(baseline = 0, contrast = 1, factor_subset = "W")
  likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "W")))
  expect_error(likelihood$get_likelihood(tmle_task, "A"))
  expect_error(likelihood$get_likelihood(tmle_task, "Y"))
})

test_that("Fitted likelihood values available for multiple factors", {
  tmle_spec <- tmle_RR(baseline = 0, contrast = 1, factor_subset = c("W", "A"))
  likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "W")))
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "A")))
  expect_error(likelihood$get_likelihood(tmle_task, "Y"))
})

test_that("Fitted likelihood values available for all factors by default", {
  tmle_spec <- tmle_RR(baseline = 0, contrast = 1)
  likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "W")))
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "A")))
  expect_true(is.numeric(likelihood$get_likelihood(tmle_task, "Y")))
})

