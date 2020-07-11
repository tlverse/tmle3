#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#

tmle3_Spec_survival <- R6Class(
  classname = "tmle3_Spec_survival",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, target_times = NULL, variable_types = NULL,...) {
      super$initialize(
        # TODO: check variable types
        # TODO: support multi-level treatments and etc
        treatment_level = treatment_level,
        control_level = control_level, 
        variable_types = variable_types, 
        target_times = target_times, 
        ...
      )
    },

    # TODO: check
    transform_data = function(data, node_list) {
      T_tilde_name <- node_list$T_tilde
      Delta_name <- node_list$Delta
      T_tilde_data <- data[T_tilde_name]
      Delta_data <- data[Delta_name]
      k_grid <- 1:max(T_tilde_data)

      if (is.null(node_list$id)) {
        id <- 1:nrow(data)
        data <- cbind(id=id, data)
        node_list$id <- "id"
      }

      all_times <- lapply(k_grid, function(t_current){
        df_time <- copy(data)
        # TODO: check
        df_time$N <- as.numeric(t_current == T_tilde_data & Delta_data == 1)
        df_time$A_c <- as.numeric(t_current == T_tilde_data & Delta_data == 0)
        df_time$pre_failure <- as.numeric(t_current <= T_tilde_data)
        df_time$t <- t_current
        return(df_time)
      })
      df_long <- rbindlist(all_times)

      long_node_list <- copy(node_list)
      long_node_list$time <- "t" 
      long_node_list$N <- "N" 
      long_node_list$A_c <- "A_c"
      long_node_list$pre_failure <- "pre_failure"

      return(list(long_data=df_long, long_node_list=long_node_list))
    },

    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      
      tmle_task <- survival_tx_task(data, node_list, survival_tx_npsem, variable_types)

      return(tmle_task)
    },

    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- survival_tx_likelihood(tmle_task, learner_list)
      return(likelihood)
    },

    make_params = function(tmle_task, likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level

      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)

      # TODO: currently support treatment specific
      # TODO: check
      param_surv <- Param_survival$new(likelihood, treatment, 
                                       target_times = self$options$target_times, 
                                       outcome_node = "N")
      tmle_params <- list(param_surv)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

# TODO
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @param target_times the time points to be targeted at during the TMLE adjustment
#' @export
# TODO: check variable types
tmle_survival <- function(treatment_level, control_level, target_times = NULL, variable_types = NULL) {
  tmle3_Spec_survival$new(treatment_level, control_level, target_times, variable_types)
}
