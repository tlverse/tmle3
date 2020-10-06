context("ATT interventions: treatment effect amongst the treated")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# setup data for test
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 1)
data$parity <- data$sexn + data$parity  * (1-data$smoked)
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0



library(simcausal)
D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = -0.8, max = 0.8) +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("A", distr = "rbinom", size = 1,  prob = plogis(W1)) +
  node("g1", distr = "rconst", const = plogis(W1)) +
  node("Y", distr = "rbinom", size =1 , prob = plogis(( -1.5 + 1 + A + W - W1/2 ))) +
  node ("EY1", distr = "rconst", const = plogis(( -1.5 + 1 + 1 + W - W1/2 ))) +
node ("EY0", distr = "rconst", const = plogis(( -1.5 + 1 + 0 + W - W1/2 )))
setD <- set.DAG(D)
data <- sim(setD, n = 10000)
data <- as.data.table(data)
print(sum((data$EY1 - data$EY0) * data$g1)/sum(data$A))


D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = -0.8, max = 0.8) +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("A", distr = "rbinom", size = 1,  prob = plogis(W1)) +
  node("g1", distr = "rconst", const = plogis(W1)) +
  node("Y", distr = "rbinom", size =1 , prob = plogis(( -1.5 + 1 + A + W - W1/2 ))) +
  node ("EY1", distr = "rconst", const = plogis(( -1.5 + 1 + 1 + W - W1/2 )))
setD <- set.DAG(D)
data <- sim(setD, n = 1000)
data <- as.data.table(data)

node_list <- list(
  W = c("W", "W1" ),
  A = "A",
  Y = "Y"
)



qlib <- make_learner(
  "Lrnr_xgboost"
)

glib <- make_learner(
  "Lrnr_xgboost"
)

logit_metalearner <- make_learner(
  Lrnr_solnp, metalearner_logistic_binomial,
  loss_loglik_binomial
)
Q_learner <- qlib
g_learner <- glib
learner_list <- list(Y = Q_learner, A = g_learner)

tmle_spec <- tmle_ATT(1, 0)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# LF_fit$undebug("get_likelihood")
# estimate likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], type = "mean")
initial_likelihood$add_factors(list(A_factor))

print(names(initial_likelihood$factor_list))

lik <- ((initial_likelihood$get_likelihood(tmle_task, "A")))
print(head(lik))
A <- tmle_task$get_tmle_node("A")
print("risk")
print(mean(-1 * ifelse(A==1, log(lik), log(1 - lik))))
print(mean(- log(lik)))
pA1 <- initial_likelihood$get_likelihoods(tmle_task, "A")
print(quantile(pA1))

updater <- tmle3_Update$new(
  cvtmle = FALSE, convergence_type = "sample_size",
  constrain_step = T, one_dimensional = T, delta_epsilon = c(-0.0005, 0.0005),
  optim_delta_epsilon = T
)
print(updater$step_number)
# debugonce(updater$update_step)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)


# define parameter
tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
att <- tmle_params[[1]]
print("before")
print(tmle_params[[1]]$estimates(tmle_task)$psi)
# fit tmle update
tmle_fit <- fit_tmle3(
  tmle_task, targeted_likelihood, list(att), updater,
  max_it
)

lik <- ((targeted_likelihood$get_likelihood(tmle_task, "A")))
A <- tmle_task$get_tmle_node("A")
print("risk")
print(mean(-1 * ifelse(A==1, log(lik), log(1 - lik))))
print(mean(- log(lik)))
print(updater$step_number)
# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se
tmle3_epsilon <- updater$epsilons[[1]]$Y

print("here")
g <-  ((targeted_likelihood$get_likelihood(tmle_task, "A")))
print(mean(tmle_params[[1]]$clever_covariates(tmle_task)$A * (g - A)))
g <- ifelse(A==1, g, 1-g)
print(mean(tmle_params[[1]]$clever_covariates(tmle_task)$A * (g - A)))
print(mean(unlist(tmle_fit$estimates[[1]]$IC)))
#################################################
# compare with the tmle package
library(tmle)

# construct likelihood estimates

# tasks for A=1 and A=0
cf_task1 <- att$cf_likelihood_treatment$cf_tasks[[1]]
cf_task0 <- att$cf_likelihood_control$cf_tasks[[1]]

# get Q
EY1 <- initial_likelihood$get_likelihoods(cf_task1, "Y")
EY0 <- initial_likelihood$get_likelihoods(cf_task0, "Y")
EY1_final <- targeted_likelihood$get_likelihoods(cf_task1, "Y")
EY0_final <- targeted_likelihood$get_likelihoods(cf_task0, "Y")
# EY0 <- rep(0, length(EY1)) # not used
Q <- cbind(EY0, EY1)

# get G
pA1 <- initial_likelihood$get_likelihoods(cf_task1, "A")
print(quantile(pA1))

# debugonce(oneStepATT)
tmle_classic_fit <- tmle(
  Y = tmle_task$get_tmle_node("Y"),
  A = tmle_task$get_tmle_node("A"),
  W = cbind(tmle_task$get_tmle_node("W"), tmle_task$get_tmle_node("W")),
  Q = Q,
  g1W = pA1,
  family = "binomial",
  alpha = 0.995,
  target.gwt = FALSE
)


# extract estimates
classic_psi <- tmle_classic_fit$estimates$ATT$psi
print( tmle3_psi)
print( classic_psi)

classic_se <- sqrt(tmle_classic_fit$estimates$ATT$var.psi)
tol <- 1 / (tmle_task$nrow)
test_that("psi matches result from classic package", {
  expect_equal(tmle3_psi, classic_psi, tol)
})

test_that("se matches result from classic package", {
  expect_equal(tmle3_se, classic_se, tol)
})
