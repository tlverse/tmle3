context("Sampling")

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

# verify we can obtain samples
samp_W_vals <- likelihood$factor_list$W$sample(tmle_task[1:50], 30, return_values=TRUE)
samp_W_task <- likelihood$factor_list$W$sample(tmle_task[1:50], 30, return_values=FALSE)

samp_A_vals <- likelihood$factor_list$A$sample(tmle_task[1:50], 30, return_values=TRUE)
samp_A_task <- likelihood$factor_list$A$sample(tmle_task[1:50], 30, return_values=FALSE)

samp_Y_vals <- likelihood$factor_list$Y$sample(tmle_task[1:50], 30, return_values=TRUE)
samp_Y_task <- likelihood$factor_list$Y$sample(tmle_task[1:50], 30, return_values=FALSE)

intervention <- define_lf(LF_static, "A", value = 1)
intervention$sample(tmle_task[1:50], 30, return_values=TRUE)
intervention$sample(tmle_task[1:50], 30, return_values=FALSE)
