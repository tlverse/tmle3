context("Caching of Likelihood objects")

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

logit_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial, loss_loglik_binomial)
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

W_factor <- define_lf(LF_emp, "W")
A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], cache = FALSE)
Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")

# construct and train likelihood
factor_list <- list(W_factor, A_factor, Y_factor)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)

# get likelihood values to populate cache
likelihood_values <- likelihood$get_likelihoods(tmle_task)

A_values <- likelihood$cache$get_values(A_factor, tmle_task, "full")
Y_values <- likelihood$cache$get_values(Y_factor, tmle_task, "full")

test_that("caching works", expect_length(Y_values, tmle_task$nrow))
test_that("disabling caching works", expect_null(A_values))
