context("Risk ratios")

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

# metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_glm_fast)
g_learner <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_RR(baseline = 0, contrast = 1)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")

targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)

# define param
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)

# extract results
summary <- tmle_fit$summary
