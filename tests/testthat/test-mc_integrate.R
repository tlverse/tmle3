context("MC integration")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)


# setup data for test
data(cpp)
data <- cpp
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0
node_list <- list(
  W = c("sexn"),
  A = "parity",
  Y = "haz01"
)

qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_xgboost"
)

logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)

mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)

Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_ATE(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# estimate likelihood
likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
# debugonce(likelihood$factor_list[["Y"]]$sample)

# TODO: rewrite sampling test below
# verify we can obtain one sample
tmle_task2 <- likelihood$sample(tmle_task$nrow)

# define param and cf likelihood
lf_static <- define_lf(LF_static, "A", value = 1)
tsm <- Param_TSM$new(likelihood, lf_static)
cf_likelihood <- tsm$cf_likelihood
fun <- function(tmle_task){mean(tmle_task$get_tmle_node("Y"))}

result <- resample(cf_likelihood, fun)

# TODO: format results in resample
result <- unlist(result)
est <- mean(result)
se <- sd(result)/sqrt(length(result))

param_ests <- tsm$estimates(tmle_task)
param_est <- param_ests$psi
ic_se <- sd(param_ests$IC)/sqrt(length(param_ests$IC))

test_that("mc integration is the same as simple mean up to mc error", 
          expect_equal(est, param_est, tol=se))
