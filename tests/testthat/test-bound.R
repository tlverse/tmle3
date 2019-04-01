context("Bounding for likelihood factors")

# generate likelihood with extreme values
set.seed(1234)
n <- 1000
W <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, plogis(10 * (W - 0.5)))
Y <- rbinom(n, 1, plogis(10 * (A - 0.5)))

data <- data.table(W, A, Y)
node_list <- list(W = "W", A = "A", Y = "Y")

tmle_spec <- tmle_TSM_all()
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

lrnr_glm <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y = lrnr_glm, A = lrnr_glm)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", learner = learner_list[["A"]]),
  define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)

# bounded_likelihood
bound_level <- 0.025
bounded_factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = bound_level),
  define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
)

bounded_likelihood_def <- Likelihood$new(bounded_factor_list)
bounded_likelihood <- bounded_likelihood_def$train(tmle_task)

g_preds <- likelihood$get_likelihood(tmle_task, "A")
bounded_preds <- bounded_likelihood$get_likelihood(tmle_task, "A")

test_that("bounds are being respected", {
  expect_lt(min(g_preds), bound_level)
  expect_gte(min(bounded_preds), bound_level)
  expect_gt(max(g_preds), 1 - bound_level)
  expect_lte(max(bounded_preds), 1 - bound_level)
})

Q_preds <- likelihood$get_likelihood(tmle_task, "Y")

updater <- tmle_spec$make_updater()
targeted_likelihood <- tmle_spec$make_targeted_likelihood(likelihood, updater)
params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- params
submodel_data <- updater$generate_submodel_data(targeted_likelihood, tmle_task)

Q_submodel <- submodel_data$Y$initial
Q_bound_level <- 1 / n
test_that("bounds are being respected in submodel", {
  expect_gt(max(Q_preds), 1 - Q_bound_level)
  expect_lte(max(Q_submodel), 1 - Q_bound_level)
})
