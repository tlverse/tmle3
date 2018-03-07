context("Shift interventions: TSM")

library(sl3)
library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)

# setup data for test
## ------------------------------------------------------------------------
# simulate simple data for tmle-shift sketch
n_obs <- 1000  # number of observations
n_w <- 1  # number of baseline covariates
tx_mult <- 2  # multiplier for the effect of W = 1 on the treatment

## baseline covariate -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, 0.5)))

## create treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))

# create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

data <- data.table(W,A,Y)
node_list <- list(W="W", A="A", Y="Y")
# setup learners
# 
# # SL learners to be used for most fits (e.g., IPCW, outcome regression)
# lrn1 <- Lrnr_mean$new()
# lrn2 <- Lrnr_glm_fast$new()
# lrn3 <- Lrnr_randomForest$new()
# sl_lrn <- Lrnr_sl$new(learners = list(lrn1, lrn2, lrn3),
#                       metalearner = Lrnr_nnls$new())
# 
# # SL learners for conditional densities to be used for the propensity score fit
# lrn1_dens <- Lrnr_condensier$new(nbins = 35, bin_estimator = lrn1,
#                                  bin_method = "equal.len")
# lrn2_dens <- Lrnr_condensier$new(nbins = 25, bin_estimator = lrn2,
#                                  bin_method = "equal.len")
# sl_lrn_dens <- Lrnr_sl$new(learners = list(lrn1_dens, lrn2_dens),
#                            metalearner = Lrnr_solnp_density$new())
# 

# use learners without randomness for ease of testing
lrn1 <- Lrnr_mean$new()

# # SL learners for conditional densities to be used for the propensity score fit
lrn1_dens <- Lrnr_condensier$new(nbins = 35, bin_estimator = lrn1,
                                 bin_method = "equal.len")

Q_learner <- lrn1
g_learner <- lrn1_dens
learner_list <- list(Y=Q_learner, A=g_learner)

tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
likelihood <- tmle_spec$make_likelihood(tmle_task, learner_list)

# define shift
# todo: put this somewhere sensible
shift_function <- function(tmle_task){
  delta=0.5
  tmle_task$get_tmle_node("A") + delta
}

shift_inverse <- function(tmle_task){
  delta=0.5
  tmle_task$get_tmle_node("A") - delta
}

# define parameter
intervention <- define_lf(LF_shift, "A", likelihood$factor_list[["A"]], shift_function, shift_inverse)
tsm <- define_param(Param_TSM, likelihood, intervention)

# define update method (submodel + loss function)
updater <- tmle_spec$make_updater(likelihood, list(tsm))

# fit tmle update
# debugonce(intervention$get_likelihood)
tmle_fit <- fit_tmle3(tmle_task, likelihood, list(tsm), updater)

# extract results
tmle3_psi <- tmle_fit$summary$tmle_est
tmle3_se <- tmle_fit$summary$se

#################################################
# compare with the txshift package
library(txshift)

# debugonce(tmle_txshift)

# todo: validate that we're getting the same errors on g fitting
tmle_sl_shift_1 <- tmle_txshift(W = W, A = A, Y = Y, delta = 0.5,
                                fluc_method = "standard",
                                ipcw_fit_args = NULL,
                                g_fit_args = list(fit_type = "sl",
                                                  sl_lrnrs = g_learner),
                                Q_fit_args = list(fit_type = "sl",
                                                  sl_lrnrs = Q_learner)
)

summary(tmle_sl_shift_1)
classic_psi <- tmle_sl_shift_1$psi
classic_se <- sqrt(tmle_sl_shift_1$var)

# only approximately equal (although it's O(1/n))
test_that("psi matches result from classic package", expect_equal(tmle3_psi, classic_psi, tol=1e-3))

# only approximately equal (although it's O(1/n))
test_that("se matches result from classic package", expect_equal(tmle3_se, classic_se, tol=1e-3))
