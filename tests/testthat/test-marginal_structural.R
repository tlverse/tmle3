context("Stratification - estimate TSM in strata")

library(sl3)
library(uuid)
library(assertthat)
library(data.table)
library(future)


# data
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0

### discrete ###
node_list <- list(
  W = c("whz"),
  V = "sexn",
  A = "parity01",
  Y = "haz01"
)

# leaners
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

# estimators
tmle_spec <- tmle_MSM()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
# disable cvtmle for this test to compare with tmle package
updater <- tmle3_Update$new(cvtmle = FALSE, convergence_type = "sample_size")
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define parameter
msm <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- msm

# fit
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, msm, updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates
cf_task1 <- msm$cf_likelihoods[["A_1"]]$cf_tasks[[1]]
cf_task0 <- msm$cf_likelihoods[["A_0"]]$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
Q <- cbind(EY0, EY1)

# get g
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")
pA0 <- initial_likelihood$get_likelihoods(cf_task0, "A")
h <- cbind(pA0, pA1)

tmle_classic_fit <- tmleMSM(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = setnames(tmle_task$get_data(NULL, node_list[["W"]]), "W"),
  V = setnames(tmle_task$get_data(NULL, node_list[["V"]]), "V"),
  MSM = "A + V",
  Q = Q,
  hAV = h,
  g1W = pA1,
  family = "binomial"
)

# extract estimates
classic_psi <- tmle_classic_fit$psi
classic_se <- tmle_classic_fit$se

# only approximately equal (although it's O(1/n))
names(tmle3_psi) <- c("A_0", "A_1", "V")
classic_psi["A"] <- classic_psi["A"] + classic_psi["(Intercept)"]
names(classic_psi) <- c("A_0", "A_1", "V")
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol = 1e-3)
})

# only approximately equal (although it's O(1/n))
tmle3_se <- tmle3_se[c(1,3)]
names(tmle3_se) <- c("A_0", "V")
classic_se <- classic_se[c(1,3)]
names(classic_se) <- c("A_0", "V")
test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol = 1e-3)
})

### continuous ###
node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  V = "agedays",
  A = "whz",
  Y = "haz"
)
processed <- process_missing(data[1:500,], node_list)
data <- processed$data
node_list <- processed$node_list

# leaners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_xgboost <- make_learner(Lrnr_xgboost)

ls_metalearner <- make_learner(Lrnr_nnls)
mn_metalearner <- make_learner(
  Lrnr_solnp, metalearner_linear_multinomial,
  loss_loglik_multinomial
)
sl_Y <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_xgboost),
  metalearner = ls_metalearner
)
sl_A <- make_learner(Lrnr_density_semiparametric)

learner_list <- list(A = sl_A, Y = sl_Y)

# estimators
tmle_spec <- tmle_MSM()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# define update method (submodel + loss function)
updater <- tmle3_Update$new()
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

# define parameter
msm <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- msm

# fit
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, msm, updater)

# extract results
tmle_ests <- tmle_fit$summary$tmle_est

