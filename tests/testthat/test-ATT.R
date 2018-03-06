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

metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, qlib, metalearner)
g_learner <- make_learner(Lrnr_sl, glib, metalearner)
learner_list <- list(Y=Q_learner, A=g_learner)
tmle_spec <- tmle_tsm_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# LF_fit$undebug("get_likelihood")
# estimate likelihood
factor_list <- list(
  define_lf(LF_np, "W", NA),
  define_lf(LF_fit, "A", type = "density", learner = learner_list[["A"]]),
  define_lf(LF_fit, "Y", type = "mean", learner = learner_list[["Y"]])
)

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(tmle_task)

# define parameter
intervention_treatment <- define_lf(LF_static, "A", value = 1)
intervention_control <- define_lf(LF_static, "A", value = 0)
att <- define_param(Param_ATT, likelihood, "Y", intervention_treatment, intervention_control)
self <- att
# define update method (submodel + loss function)
updater <- tmle_spec$make_updater(likelihood, list(tsm))

# fit tmle update
tmle_fit <- fit_tmle3(tmle_task, likelihood, list(tsm), updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# task for A=1
cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))

# get Q
EY1 <- likelihood$get_initial_likelihoods(cf_task, "Y")
EY1_final <- likelihood$get_likelihoods(cf_task, "Y")
EY0 <- rep(0, length(EY1)) # not used
Q = cbind(EY0, EY1)

# get G
pA1 <- likelihood$get_initial_likelihoods(cf_task, "A")
pDelta1 <- cbind(pA1, pA1)
tmle_classic_fit <- tmle(Y=tmle_task$get_tmle_node("Y"), 
                            A=NULL,
                            W=tmle_task$get_tmle_node("W"),
                            Delta = tmle_task$get_tmle_node("A"),
                            Q = Q,
                            pDelta1 = pDelta1)

# extract estimates
classic_psi <- tmle_classic_fit$estimates$EY1$psi
classic_se <- sqrt(tmle_classic_fit$estimates$EY1$var.psi)

# only approximately equal (although it's O(1/n))
test_that("psi matches result from classic package", expect_equal(tmle3_psi, classic_psi, tol=1e-3))

# only approximately equal (although it's O(1/n))
test_that("se matches result from classic package", expect_equal(tmle3_se, classic_se, tol=1e-3))
