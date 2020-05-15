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
    initialize = function(treatment_level, control_level, ...) {
      super$initialize(
        # TODO: support multi-level treatments and etc
        treatment_level = treatment_level,
        control_level = control_level, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types

      # TODO: random initialize dN and dA_c as binomials; function design
      # data[, "dN"] <- data[[node_list$T_tilde]]
      # data[, "dA_c"] <- data[[node_list$T_tilde]]

      # data[, "dN"] <- rbinom(nrow(data), 1,.5)
      # data[, "dA_c"] <- rbinom(nrow(data), 1,.5)
      # node_list["dN"] <- "dN"
      # node_list["dA_c"] <- "dA_c"
      
      tmle_task <- survival_tx_task(data, node_list, survival_tx_npsem, variable_types)

      return(tmle_task)
    },
    # TODO
    # make_initial_likelihood = function(tmle_task, learner_list = NULL) {
    #   # produce trained likelihood when likelihood_def provided

    #   if (!is.null(self$options$likelihood_override)) {
    #     likelihood <- self$options$likelihood_override$train(tmle_task)
    #   } else {
    #     likelihood <- survival_tx_likelihood(tmle_task, learner_list)
    #   }

    #   return(likelihood)
    # }

    # make_long_tmle_task = function(long_data, long_node_list, ...) {
    #   variable_types <- self$options$variable_types

    #   # TODO: function design
    #   long_tmle_task <- survival_tx_task(long_data, long_node_list, survival_tx_long_npsem, variable_types)

    #   return(long_tmle_task)
    # },

    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided

      # TODO: check if needed
      # if (!is.null(self$options$likelihood_override)) {
      #   base_likelihood <- self$options$likelihood_override$train(tmle_task)
      # } else {
      #   # make likelihood only for W, A
      #   base_likelihood <- survival_tx_base_likelihood(tmle_task, learner_list)
      # }

      # # make likelihood only for W, A
      # base_likelihood <- survival_tx_base_likelihood(tmle_task, learner_list)

      # # generate likelihood estimates for W, A
      # base_likelihood$get_likelihoods(tmle_task)

      # # transform original data into long version
      # short_data <-tmle_task$data
      # short_npsem <- tmle_task$npsem
      # long_data <- make_long_data(short_data, short_npsem)
      # long_node_list <- make_long_node_list(short_npsem)

      # # make long tmle task
      # selected_long_data <- long_data[which(long_data["N"] == 0 & long_data["A_c"] == 0),]
      # long_tmle_task <- self$make_long_tmle_task(selected_long_data, long_node_list)

      # # make likelihood for dN and dA_c
      # likelihood <- survival_tx_hazard_likelihood(long_tmle_task, learner_list)

      # # generate likelihood estimates for dN and dA_c
      # full_tmle_task <- self$make_long_tmle_task(long_data, long_node_list)
      # # TODO: set rs[72] <- 1 - rs[72] 
      # likelihood$get_likelihoods(full_tmle_task)

      # TODO: merge likelihoods
      likelihood <- survival_tx_likelihood(tmle_task, learner_list)
      return(likelihood)
    }
  ),
  active = list(),
  private = list()
)

# TODO
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_survival <- function(treatment_level, control_level) {
  tmle3_Spec_survival$new(treatment_level, control_level)
}
