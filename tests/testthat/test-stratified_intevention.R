context("Stratification - estimate TSM in strata")

library(sl3)
library(uuid)
library(assertthat)
library(data.table)
library(future)

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
  V = "sexn",
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

logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
ate_spec <- tmle_ATE(1, 0)
strat_spec <- tmle_stratified(ate_spec, "mrace")
tmle_spec <- strat_spec

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = TRUE)

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

params <- tmle_spec$make_params(tmle_task, targeted_likelihood)

tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, params, updater)
tmle_fit <- tmle3(strat_spec, data, node_list, learner_list)
