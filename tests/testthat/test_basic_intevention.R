library(tmle3)
library(testthat)
library(sl3)
library(uuid)
library(assertthat)
library(data.table)
context("test_basic_intervention.R -- TSM and ATE for single node interventions")


data(cpp)
cpp <- cpp[!is.na(cpp[, "haz"]), ]
cpp$parity3 <- cpp$parity
cpp$parity3[cpp$parity3>3] <- 3
cpp$parity01 <- as.numeric(cpp$parity>0)
cpp[is.na(cpp)] <- 0
cpp$haz01 <- as.numeric(cpp$haz>0)
tmle_nodes <- list(define_node("W", c("apgar1", "apgar5", "gagebrth", "mage", "meducyrs", "sexn")),
                   define_node("A", c("parity01"), c("W")),
                   define_node("Y", c("haz01"), c("A","W")))

task <- tmle_Task$new(cpp, tmle_nodes=tmle_nodes)

glm_fast <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)

# define and fit likelihood
factor_list <- list(define_lf(LF_static,"W", NA),
                    define_lf(LF_fit,"A", glm_fast),
                    define_lf(LF_fit,"Y", lrnr_mean))

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(task)

# define parameter and get tmle likelihood
intervention = define_cf(define_lf(LF_static, "A", value=1))

tsm <- Param_TSM$new(intervention)
tmle_likelihood <- fit_tmle_likelihood(likelihood, task, tsm)
init_ests <- tsm$estimates(likelihood, task)
tmle_ests <- tsm$estimates(tmle_likelihood, task)

# get initial and parameter estimates
tmle3_init_psi <- init_ests$psi
tmle3_tmle_psi <- tmle_ests$psi
ED <- mean(tmle_ests$IC)
expect_lt(abs(ED), 1/task$nrow)


# compare to gentmle2
library(gentmle2)
Q1W = likelihood$get_predictions(cf_task_1, "Y")
QAW = likelihood$get_predictions(task, "Y")
g1W = likelihood$get_predictions(task, "A")
A <- as.matrix(task$get_regression_task("A")$Y)
Y <- as.matrix(task$get_regression_task("Y")$Y)
# this tmle doesn't need Q0k but because it's hard-coded in gentmle we need to provide it
tmledata <- data.frame(A=A, Y=Y, gk=g1W, Qk = QAW, Q1k=Q1W, Q0k=Q1W)

result <- gentmle(tmledata, params=list(EY1=param_EY1), 
                  submodel = Q_submodel_logit, loss = Q_loss_loglik)
gentmle2_init_est <- result$initests
gentmle2_tmle_est <- result$tmleests

expect_equivalent(tmle3_init_psi, gentmle2_init_est)
expect_equal(tmle3_tmle_psi, gentmle2_tmle_est, tolerance = 1e-4, check.attributes = FALSE)


cf_a1 = define_cf(define_lf(LF_static, "A", value=1))
cf_a0 = define_cf(define_lf(LF_static, "A", value=0))
ate <- Param_ATE$new(cf_a0, cf_a1)

tmle_likelihood <- fit_tmle_likelihood(likelihood, task, ate)

init_ests <- ate$estimates(likelihood, task)
init_ests$psi
tmle_ests <- ate$estimates(tmle_likelihood, task)
tmle_ests$psi
ED <- mean(tmle_ests$IC)
expect_lt(abs(ED), 1/task$nrow)
