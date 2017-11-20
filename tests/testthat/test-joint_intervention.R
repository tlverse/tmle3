context("Joint interventions: TSM and ATE for joint interventions")

library(sl3)
library(uuid)
library(assertthat)
library(data.table)

data(cpp)
cpp <- cpp[!is.na(cpp[, "haz"]), ]
cpp$parity3 <- cpp$parity
cpp$parity3[cpp$parity3>3] <- 3
cpp$parity01 <- as.numeric(cpp$parity>0)
cpp[is.na(cpp)] <- 0
cpp$haz01 <- as.numeric(cpp$haz>0)
cpp$sex01 <- as.numeric(cpp$sexn-1)
tmle_nodes <- list(define_node("W", c("apgar1", "apgar5", "gagebrth", "mage",
                                      "meducyrs", "sexn")),
                   define_node("A", c("parity01"), c("W")),
                   define_node("Z", c("sex01"), c("A", "W")),
                   define_node("Y", c("haz01"), c("A","A2","W")))

task <- tmle_Task$new(cpp, tmle_nodes=tmle_nodes)
glm_fast <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)

# define and fit likelihood
factor_list <- list(define_lf(LF_static,"W", NA),
                    define_lf(LF_fit,"A", glm_fast),
                    define_lf(LF_fit,"Z", glm_fast),
                    define_lf(LF_fit,"Y", glm_fast))

likelihood_def <- Likelihood$new(factor_list)
likelihood <- likelihood_def$train(task)

# define parameter and get tmle likelihood
intervention = define_cf(c(define_lf(LF_static, "A", value = 1),
                           define_lf(LF_static, "Z", value = 1)))

tsm <- Param_TSM$new(intervention)
tmle_likelihood <- fit_tmle_likelihood(likelihood, task, tsm)

init_ests <- tsm$estimates(likelihood, task)
init_ests$psi
tmle_ests <- tsm$estimates(tmle_likelihood, task)
tmle_ests$psi

# TEST: mean of EIF is nearly zero.
ED <- mean(tmle_ests$IC)
expect_lt(abs(ED), 1 / task$nrow)

# controlled direct effect
a0z0 = define_cf(c(define_lf(LF_static, "A", value = 0),
                   define_lf(LF_static, "Z", value = 0)))

a1z0 = define_cf(c(define_lf(LF_static, "A", value = 1),
                   define_lf(LF_static, "Z", value = 0)))

cde <- Param_ATE$new(a0z0, a1z0)

tmle_likelihood <- fit_tmle_likelihood(likelihood, task, cde)

init_ests <- cde$estimates(likelihood, task)
init_ests$psi
tmle_ests <- cde$estimates(tmle_likelihood, task)
tmle_ests$psi

# TEST: mean of EIF is nearly zero.
ED <- mean(tmle_ests$IC)
expect_lt(abs(ED), 1 / task$nrow)

