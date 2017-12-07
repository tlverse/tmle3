#' Class for Estimation of ATE Parameter
#'
#' A class for estimating the Average Treatment Effect (ATE) parameter.
#'
#' @importFrom R6 R6Class
#'
#' @export
#
Param_ATE <- R6Class(classname = "Param_ATE",
                     portable = TRUE,
                     class = TRUE,
                     inherit = Param_base,
  public = list(
    initialize = function(counterfactual_0, counterfactual_1,
                           outcome_node = "Y") {
      private$.counterfactual_0 = counterfactual_0
      private$.counterfactual_1 = counterfactual_1
      private$.outcome_node = outcome_node
    },
    HA = function(likelihood, task) {
      tsm_1 = Param_TSM$new(self$counterfactual_1, self$outcome_node)
      tsm_0 = Param_TSM$new(self$counterfactual_0, self$outcome_node)
      HA = tsm_1$HA(likelihood, task) - tsm_0$HA(likelihood, task)
      return(HA)
    },
    estimates = function(likelihood, task) {
      cf_task_1 = self$counterfactual_1$cf_task(task)
      cf_task_0 = self$counterfactual_0$cf_task(task)

      Y = task$get_regression_task(self$outcome_node)$Y
      HA = self$HA(likelihood, task)

      EY = likelihood$get_predictions(task, self$outcome_node)
      EY1 = likelihood$get_predictions(cf_task_1, self$outcome_node)
      EY0 = likelihood$get_predictions(cf_task_0, self$outcome_node)

      psi = mean(EY1 - EY0)
      IC = HA * (Y - EY) + EY1 - EY0 - psi

      result = list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    counterfactual_0 = function() {
      return(private$.counterfactual_0)
    },
    counterfactual_1 = function() {
      return(private$.counterfactual_1)
    }
  ),
  private = list(
    .counterfactual_0 = NULL,
    .counterfactual_1 = NULL
  )
)

