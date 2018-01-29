context("Basic interventions: TSM and ATE for single node interventions.")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)


# setup data for test
data(cpp)
data <- as.data.table(cpp)
# data$parity01 <- as.numeric(data$parity > 0)
# data$parity01_fac <- factor(data$parity01)
# data$haz01 <- as.numeric(data$haz > 0)

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

#todo: fix whater's going on here with binomial
#some failure with predictions I think
#also failure of get_levels
#probably want to use values insteads
#todo: figure out how to correct for Ystar
# lrnr_mean <- make_learner("Lrnr_mean")
qlib <- glib <- make_learner_stack("Lrnr_glmnet")
mn_metalearner <- make_learner(Lrnr_solnp, loss_function = loss_loglik_multinomial, learner_function = metalearner_linear_multinomial)
metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, qlib, metalearner)
g_learner <- make_learner(Lrnr_sl, glib, mn_metalearner)
learner_list <- list(Y=Q_learner, A=g_learner)
tmle_fit <- tmle3(tmle_tsm_all(), data, node_list, learner_list)
tmle_fit$summary
tmle3_Spec$debug("make_likelihood")
tmle3_Update$debug("update_step")




psis <- tmle_fit$psi
ICs <- tmle_fit$IC
psi_baseline <- psis[baseline]
psi_other <- psis[-baseline]
IC_baseline <- ICs[,baseline]
IC_other <- ICs[,-baseline]
psi_contrasts <- psi_other - psi_baseline
IC_contrasts <- apply(IC_other, 2, function(IC) IC-IC_baseline)

# debugonce(tmle_param$estimates)
# debugonce(tmle_param$cf_likelihood$get_likelihoods)
tmle_task <- tmle_fit$tmle_task
test <- tmle_param$estimates(tmle_task)
plot(tmle_fit)

# unbound <- function(x, lower, upper){
#   (upper - lower) * x + lower
# }
# 
# unbound(tmle_fit$psi, lower, upper)
