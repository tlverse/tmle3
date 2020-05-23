library(tmle3)
library(sl3)
library(data.table)
library(R6)
library(SuperLearner)
library(origami)
library(uuid)

vet_data <- read.csv("https://raw.githubusercontent.com/tlverse/deming2019-workshop/master/data/veteran.csv")
vet_data$trt <- vet_data$trt - 1

# TODO: check
vet_data$time <- ceiling(vet_data$time/20) # make fewer times for testing
tmax <- max(vet_data$time)
all_times <- lapply(seq_len(tmax), function(t_current){
  vet_data_time <- copy(vet_data)
  # TODO: check
  vet_data_time$N <- ifelse(t_current == vet_data_time$time & vet_data_time$status == 1, 1, 0)
  vet_data_time$A_c <- ifelse(t_current == vet_data_time$time & vet_data_time$status == 0, 1, 0)
  vet_data_time$t <- t_current

  return(vet_data_time)
})
vet_data_long <- rbindlist(all_times)

# N <- ifelse(vet_data$time <= vet_data$t & vet_data$status == 1, 1, 0)
# vet_data$dN <- ifelse(N == 1 & vet_data$N == 0, 1, 0)
# A_c <- ifelse(vet_data$time <= vet_data$t & vet_data$status == 0, 1, 0)
# vet_data$dA_c <- ifelse(A_c == 1 & vet_data$A_c == 0, 1, 0)

node_list <- list(W = c("celltype", "karno", "diagtime", "age", "prior"), A = "trt", T_tilde = "time", Delta = "status", 
	t = "t", N = "N", A_c = "A_c")

lrnr_xgb <- make_learner(Lrnr_xgboost)
learner_list <- list(A = lrnr_xgb, W = lrnr_xgb, N = lrnr_xgb, A_c = lrnr_xgb)

survival_spec <- tmle_survival(treatment_level = 1, control_level = 0)
survival_task <- survival_spec$make_tmle_task(vet_data_long, node_list)

likelihood <- survival_spec$make_initial_likelihood(survival_task, learner_list)

# Update Process
initial_likelihood <- likelihood
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
tmle_task <- survival_task
tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)
tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
)

# # get likelihoods for A, W, N, A_c
# likelihood$get_likelihoods(survival_task)

# intervention_list <- define_lf(LF_static, "A", value = 1)
# observed_likelihood <- likelihood
# cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
# tmle_task <- survival_task
# fold_number <- "full"
# intervention_nodes <- "A"
# survival_param <- Param_survival$new(observed_likelihood, intervention_list, outcome_node = "N")

# # get likelihoods for the conditional hazards 
# # dN, dA_c represents hazards of failure and censoring respectively, 
# # and dN_1 is the conditional hazards at time1
# likelihood$get_hazards(survival_task)

# # get likelihoods for the survival 
# # S_N, S_A_c represents survivals of failure and censoring respectively
# likelihood$get_survival(survival_task)
