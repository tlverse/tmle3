#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_coxph <- R6Class(
  classname = "tmle3_Spec_coxph",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(formula = ~1, treatment_level, control_level, variable_types = NULL, delta_epsilon = 0.05, ...) {
      super$initialize(
        formula = formula,
        treatment_level = treatment_level,
        control_level = control_level,
        delta_epsilon = delta_epsilon,
        variable_types = variable_types,

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
        data <- cbind(id = id, data)
        node_list$id <- "id"
      }

      all_times <- lapply(k_grid, function(t_current) {
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

      return(list(long_data = df_long, long_node_list = long_node_list))
    },

    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      data_list <- self$transform_data(data, node_list)
      tmle_task <- survival_tx_task(data_list$long_data, data_list$long_node_list, variable_types)

      return(tmle_task)
    },

    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- survival_tx_likelihood(tmle_task, learner_list)
      return(likelihood)
    },
    make_updater = function(convergence_type = "sample_size", verbose = TRUE, ...) {
      if (!is.null(self$options$verbose)) {
        verbose <- self$options$verbose
      }

      updater <- tmle3_Update$new(maxit = 100, one_dimensional = TRUE, convergence_type = convergence_type, verbose = verbose, delta_epsilon = self$options$delta_epsilon, constrain_step = TRUE, bounds = c(0.0025), ...)

      return(updater)
    },

    make_params = function(tmle_task, likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level

      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)

      # TODO: currently support treatment specific
      # TODO: check
      param_surv <- Param_coxph$new(likelihood, self$options$formula, treatment, control,
        outcome_node = "N"
      )
      tmle_params <- list(param_surv)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)
