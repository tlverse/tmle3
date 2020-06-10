library(here)
source(here("Moss-simulation-master/code_simulation/simulate_data.R"))

devtools::load_all("Moss")
devtools::load_all("tmle3")
library(sl3)
library(ggplot2)

reshape_long_data = function(long_data, t_max) {
      n <- length(long_data) / t_max
      # TODO: assume long_data is a list
      rs <- list()
      for (i in 1:t_max) {
        current <- long_data[seq(1 + (i - 1) * n, i * n)]
        rs <- c(rs, list(current))
      }
      rs <- do.call(cbind, rs)
      return(rs)
}

l2_diff = function(x, y) {
  return(sqrt(sum((x - y)^2)) / length(x))
}

# context("methods for the survival curve")
# set.seed(11)
# source("./simulate_data.R")
# simulation

# TODO: set seed
set.seed(11)
n_sim <- 1e2
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat
true_surv <- simulated$true_surv1

# TODO: check
while (all(df$Delta == 1)) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1
}

# TODO: check
# set.seed(42)
# n_sim <- 1e2
# for (i in 1:5) {
#   simulated <- simulate_data(n_sim = n_sim)
#   df <- simulated$dat
#   print(all(df$Delta == 1))
# }

# TODO: lrnr
# sl_lib_g <- c("SL.xgboost")
# sl_lib_censor <- c("SL.xgboost")
# sl_lib_failure <- c("SL.xgboost")
# sl_lib_g <- c("SL.glm")
# sl_lib_censor <- c("SL.glm")
# sl_lib_failure <- c("SL.glm")
sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam")
range(df$T.tilde)
df$T.tilde <- df$T.tilde + 1
k_grid <- 1:max(df$T.tilde)

sl_fit <- initial_sl_fit(
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  A = df$A,
  W = data.frame(df[, c("W", "W1")]),
  t_max = max(df$T.tilde),
  sl_treatment = sl_lib_g,
  sl_censoring = sl_lib_censor,
  sl_failure = sl_lib_failure
)
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

test_that("sl_1 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_1$survival, is.na)))
})
test_that("sl_0 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_0$survival, is.na)))
})

################################################################################
# tlverse
tmax <- max(df$T.tilde)
all_times <- lapply(seq_len(tmax), function(t_current){
  df_time <- copy(df)
  # TODO: check
  df_time$N <- ifelse(t_current == df$T.tilde & df$Delta == 1, 1, 0)
  df_time$A_c <- ifelse(t_current == df$T.tilde & df$Delta == 0, 1, 0)
  df_time$t <- t_current

  return(df_time)
})
df_long <- rbindlist(all_times)

node_list <- list(W = c("W", "W1"), A = "A", T_tilde = "T.tilde", Delta = "Delta", 
  t = "t", N = "N", A_c = "A_c")

# TODO: check
# lrnr_xgb <- make_learner(Lrnr_xgboost)
# learner_list <- list(A = lrnr_xgb, N = lrnr_xgb, A_c = lrnr_xgb)
# lrnr_glm <- make_learner(Lrnr_glm)
# learner_list <- list(A = lrnr_glm, N = lrnr_glm, A_c = lrnr_glm)
# lrnr_mean <- make_learner(Lrnr_mean)
# learner_list <- list(A = lrnr_mean, N = lrnr_mean, A_c = lrnr_mean)
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
# lrnr_earth <- make_learner(Lrnr_earth)
# sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam, lrnr_earth))
sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(A = sl_A, N = sl_A, A_c = sl_A)

# TODO: check
var_types <- list(T_tilde = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), 
  Delta = Variable_Type$new("binomial"))
survival_spec <- tmle_survival(treatment_level = 1, control_level = 0, variable_types = var_types)
survival_task <- survival_spec$make_tmle_task(df_long, node_list)

likelihood <- survival_spec$make_initial_likelihood(survival_task, learner_list)

initial_likelihood <- likelihood
# TODO: check
# up <- tmle3_Update_survival$new(maxit = 1e2, clipping = 1e-1 / sqrt(n_sim))
# up <- tmle3_Update_survival$new(
#     one_dimensional = TRUE, constrain_step = TRUE,
#     maxit = 1e2, cvtmle = TRUE,
#     convergence_type = "sample_size",
#     delta_epsilon = 1e-2,
#     fit_method = "classic"
#   )
up <- tmle3_Update_survival$new(
  # TODO: check
  # one_dimensional = TRUE, constrain_step = TRUE,
  maxit = 3e1, cvtmle = TRUE,
  convergence_type = "sample_size",
  # delta_epsilon = 1e-2,
  fit_method = "l2",
  clipping = 1e-2
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
tmle_task <- survival_task
tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

ps <- tmle_params[[1]]
cf_task <- ps$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
pN1 <- ps$observed_likelihood$get_likelihoods(cf_task, "N")
pS_N1 <- ps$hazards_to_survival(pN1, tmax)
r_pN1 <- reshape_long_data(pN1, tmax)
r_pS_N1 <- reshape_long_data(pS_N1, tmax)
# sum(sl_fit$density_failure_1$survival - r_pS_N1)

# TODO: check
t_surv1 <- simulated$true_surv1(k_grid - 1)
# survival_truth_1 <- survival_curve$new(t = k_grid, survival = t_surv1)

psi0_moss <- colMeans(sl_fit$density_failure_1$survival)
psi0_tl <- ps$get_psi(pS_N1, tmax)
l2_diff(psi0_moss, t_surv1)
l2_diff(psi0_tl, t_surv1)

dt <- data.table(psi0_moss, psi0_tl, t_surv1)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()

################################################################################
# moss hazard submodel
moss_hazard_l2 <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
# TODO: check
rs_moss <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-1 / sqrt(n_sim), max_num_interation = 1e2, verbose = FALSE
)
psi1_moss <- rs_moss$psi_n
eic_moss <- rs_moss$eic_list

# tlverse update process
tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
)
rs <- tmle_fit_manual$estimates[[1]]
psi1_tl <- rs$psi
# sum(psi1_moss - psi1_tl)

# TODO: check
eic_tl <- up$update(targeted_likelihood, tmle_task)

l2_diff(psi1_moss, t_surv1)
l2_diff(psi1_tl, t_surv1)

dt <- data.table(psi1_moss, psi1_tl, t_surv1)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()
