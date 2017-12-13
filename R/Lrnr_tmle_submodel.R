#' Fitting submodel
#'
#' @importFrom R6 R6Class
#'
#' @export
#
Lrnr_tmle_update <- R6Class(
  classname = "Lrnr_tmle_update", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(tmle_param, lrnr_submodel, ...) {
      # this captures all parameters to initialize and saves them as self$params
      params <- args_to_list()
      super$initialize(params = params, ...)
    },

    prepare_task = function(likelihood, task) {
      # todo: modify this to support fluctuations where covariate is used as weights
      param <- self$params$tmle_param
      outcome_node <- param$outcome_node
      Y <- task$get_tmle_node(outcome_node)
      EY <- likelihood$get_predictions(task, outcome_node)
      HA <- param$HA(likelihood, task)

      tmle_data <- data.table(HA = HA, EY = EY, Y = Y)
      fluc_task <- sl3_Task$new(
        tmle_data, outcome = "Y", offset = "EY",
        covariates = "HA"
      )
      return(fluc_task)
    }
  ),
  private = list(
    .properties = c(""),

    .train = function(task) {
      tmle_fluc <- self$params$lrnr_submodel
      return(tmle_fluc$train(task))
    },

    .predict = function(task = NULL) {
      return(self$fit_object$predict(task))
    }
  )
)
