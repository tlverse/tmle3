context("Basic interventions: TSM and ATE for single node interventions.")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)


# setup data for test
data(cpp)
data <- as.data.table(cpp)

# todo: handle missingness better
data[is.na(data)] <- 0
node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  A = "parity",
  Y = "haz"
)
A_node <- node_list$A
A_vals <- unlist(data[, A_node, with = FALSE])

# todo: don't mutate cols, add new cols
if (!is.factor(A_vals)) {

  # todo: make cuts script_param
  quants <- seq(from = 0, to = 1, length = 5)
  cuts <- quantile(A_vals, quants)
  A_factor <- cut(A_vals, cuts, right = FALSE, include.lowest = TRUE)
  data[, A_node] <- A_factor
}

Y_node <- node_list$Y
Y_vals <- unlist(data[, Y_node, with = FALSE])
min_Y <- min(Y_vals)
max_Y <- max(Y_vals)
range <- max_Y - min_Y
lower <- min_Y - 0.1 * range
upper <- max_Y + 0.1 * range

Y_bounded <- (Y_vals - lower) / (upper - lower)
data[, Y_node] <- Y_bounded
qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glmnet",
  "Lrnr_glm_fast",
  "Lrnr_xgboost",
  "Lrnr_randomForest"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glmnet",
  "Lrnr_xgboost"
)

#todo: figure out how to correct for Ystar
qlib <- glib <- make_learner_stack("Lrnr_mean")
mn_metalearner <- make_learner(Lrnr_solnp, loss_function = loss_loglik_multinomial, learner_function = metalearner_linear_multinomial)
metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, qlib, metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
learner_list <- list(Y=Q_learner, A=g_learner)
tmle_fit <- tmle3(tmle_tsm_all(), data, node_list, learner_list)
tmle_fit$summary
tmle_fit$tmle_params[[1]]$cf_likelihood$name
self <- tmle_fit$tmle_params[[1]]$cf_likelihood

tmle_fit$psi
tmle_fit$ci
unbound <- function(x, lower, upper){
  (upper - lower) * x + lower
}

unbound(tmle_fit$psi, lower, upper)
