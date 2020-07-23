library(tmle3)
library(sl3)
library(data.table)

vet_data <- read.csv("https://raw.githubusercontent.com/tlverse/deming2019-workshop/master/data/veteran.csv")
vet_data$trt <- vet_data$trt - 1
# make fewer times for illustration
vet_data$time <- ceiling(vet_data$time / 20) 
k_grid <- 1:max(vet_data$time)

all_times <- lapply(k_grid, function(t_current){
  df_time <- copy(vet_data)
  # TODO: check
  df_time$N <- as.numeric(t_current == vet_data$time & vet_data$status == 1)
  df_time$A_c <- as.numeric(t_current == vet_data$time & vet_data$status == 0)
  df_time$pre_failure <- as.numeric(t_current <= vet_data$time)
  df_time$t <- t_current
    
  return(df_time)
})

df_long <- rbindlist(all_times)

node_list <- list(W = c("celltype", "karno", "diagtime", "age", "prior"), A = "trt", T_tilde = "time", Delta = "status", 
	time = "t", N = "N", A_c = "A_c", id ="X", pre_failure = "pre_failure")

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_gam <- make_learner(Lrnr_gam)
sl_A <- Lrnr_sl$new(learners = list(lrnr_mean, lrnr_glm, lrnr_gam))
learner_list <- list(A = sl_A, N = sl_A, A_c = sl_A)

var_types <- list(T_tilde = Variable_Type$new("continuous"), t = Variable_Type$new("continuous"), Delta = Variable_Type$new("binomial"))
survival_spec <- tmle_survival(treatment_level = 1, control_level = 0,
                                 variable_types = var_types)
tmle_task <- survival_spec$make_tmle_task(df_long, node_list)
initial_likelihood <- survival_spec$make_initial_likelihood(tmle_task, learner_list)

up <- tmle3_Update_survival$new(
    maxit = 3e1,
    cvtmle = TRUE,
    convergence_type = "scaled_var",
    delta_epsilon = 1e-2,
    fit_method = "l2",
    use_best = TRUE,
    verbose=TRUE
  )

targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
tmle_params <- survival_spec$make_params(tmle_task, targeted_likelihood)

max(abs(colMeans(tmle_params[[1]]$estimates(tmle_task, "validation")$IC[,1:10])))

# debugonce(tmle_params[[1]]$estimates)
tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)

# conv <- apply(abs(do.call(rbind,up$EDs)),1,max)
tmle_fit_manual$estimates[[1]]$psi

print(tmle_fit_manual)

