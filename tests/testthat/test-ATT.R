context("Basic interventions: TSM and ATE for single node interventions.")

library(sl3)
library(tmle3)
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

# LF_fit$undebug("get_likelihood")
# estimate likelihood
factor_list <- list(
  define_lf(LF_emp, "W", NA),
  define_lf(LF_fit, "A", type = "density", learner = learner_list[["A"]]),
  define_lf(LF_fit, "Y", type = "mean", learner = learner_list[["Y"]])
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)

updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

# define parameter
intervention_treatment <- define_lf(LF_static, "A", value = 1)
intervention_control <- define_lf(LF_static, "A", value = 0)
att <- define_param(
  Param_ATT, targeted_likelihood, intervention_treatment,
  intervention_control
)
updater$tmle_params <- list(att)


# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(att), updater)

# extract results
# debugonce(summary_from_estimates)
tmle_fit$summary

tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
