context("CV-TMLE: TSM for single static intervention.")

library(sl3)
# library(tmle3)
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
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)


EY <- initial_likelihood$get_likelihoods(tmle_task, "Y", 0)
EY_direct <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, 0)
all.equal(EY,EY_direct)

# define parameter
# cf_likelihood <- CF_Likelihood$new(likelihood, intervention)

# define update method (submodel + loss function)
# updater <- tmle_spec$make_updater(likelihood, list(tsm))
updater <- tmle3_cv_Update$new()

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)

intervention <- define_lf(LF_static, "A", value = 1)

# todo: make params not store likelihood info internally!
tsm <- define_param(Param_TSM, targeted_likelihood, intervention)
updater$tmle_params <- tsm



# debug(targeted_likelihood$get_likelihood)
mean(tsm$estimates(tmle_task)$psi)
# debug(targeted_likelihood$update)

# debugonce(targeted_likelihood$update)
tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, list(tsm), updater)
mean(targeted_likelihood$get_likelihoods(tmle_task, "Y"))

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# task for A=1
# cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A = 1))
cf_task <- tsm$cf_likelihood$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task, "Y", 0)
EY <- initial_likelihood$get_likelihoods(tmle_task, "Y", 0)
EY1_direct <- initial_likelihood$factor_list[["Y"]]$get_likelihood(cf_task, 0)
EY_direct <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, 0)

EY1_final <- targeted_likelihood$get_likelihoods(cf_task, "Y", 0)
EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

submodel_data <- updater$generate_submodel_data(initial_likelihood, tmle_task, 0)
z=submodel_data$Y$observed*submodel_data$Y$initial
z2=tmle_task$get_tmle_node("Y")*EY1
all.equal(z,z2)
head(cbind(z,z2))

EY <- initial_likelihood$get_likelihoods(tmle_task, "Y", 0)

all.equal(submodel_data$Y$observed, tmle_task$get_tmle_node("Y"))
all.equal(submodel_data$Y$initial, EY)
A1=which(tmle_task$get_tmle_node("A")==1)
table(match=EY==EY1,A1=tmle_task$get_tmle_node("A")==1)
initial_likelihood$cache$cache
(EY==EY1)[Y1]
length(EY1)
# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task, "A", 0)
pDelta1 <- cbind(pA1, pA1)
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = NULL,
  W = tmle_task$get_tmle_node("W"),
  Delta = tmle_task$get_tmle_node("A"),
  Q = Q,
  pDelta1 = pDelta1
)

tmle_classic_fit$epsilon
updater$epsilons
