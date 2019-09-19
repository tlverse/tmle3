context("Variable importance measure framework")

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
    "meducyrs"
  ),
  A = c("sexn", "parity01", "mracen", "smoked"),
  Y = "haz01"
)

q_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

g_lib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, q_lib, metalearner)
g_learner <- make_learner(Lrnr_sl, g_lib, metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_PAR(baseline_level = 1)

vim_results <- tmle3_vim(tmle_spec, data, node_list, learner_list)

# plot_vim(vim_results)
