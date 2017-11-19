#' Class for Estimation of TSM Parameter
#'
#' @importFrom R6 R6Class
#'
#' @export
#
Param_TSM <- R6Class(classname = "Param_TSM",
                      portable = TRUE,
                      class = TRUE,
                      inherit = Param_base,
  public = list(
    initialize <- function(counterfactual, outcome_node = "Y") {
      private$.counterfactual = counterfactual
      private$.outcome_node = outcome_node
    },
    HA <- function(likelihood, task) {
      intervention_nodes <- self$counterfactual$intervention_nodes
      cf_likelihood <- self$counterfactual$cf_likelihood(likelihood)

      pA <- likelihood$get_likelihoods(task, intervention_nodes)
      cf_pA <- cf_likelihood$get_likelihoods(task, intervention_nodes)
      HA <- cf_pA/pA

      # collapse across multiple intervention nodes
      if (!is.null(ncol(HA)) && ncol(HA) > 1) {
        HA <- apply(HA, 1, prod)
      }
      return(HA)
    },
    estimates <- function(likelihood, task) {
      cf_task <- self$counterfactual$cf_task(task)

      Y <- task$get_regression_task(self$outcome_node)$Y
      HA <- self$HA(likelihood, task)

      EY <- likelihood$get_predictions(task, self$outcome_node)
      EY1 <- likelihood$get_predictions(cf_task, self$outcome_node)

      psi <- mean(EY1)
      IC = HA * (Y - EY) + EY1 - psi

      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    counterfactual <- function() {
      return(private$.counterfactual)
    }
  ),
  private = list(
    .counterfactual <- NULL
  )
)

