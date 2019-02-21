#' Mean of Outcome Node
#'
#' Parameter for marginal mean of Y: \eqn{\Psi=E[Y]}. No TMLE update needed, but can be used in delta method calculations.
#' Useful for example, in calculating attributable risks.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_TSM, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_mean <- R6Class(
  classname = "Param_mean",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, ..., outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      return(list(Y = rep(1, tmle_task$nrow)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)

      Y <- tmle_task$get_tmle_node(self$outcome_node)
      # todo: separate out psi
      # todo: make this a function of f(W)
      psi <- mean(EY)
      IC <- Y - psi


      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[%s]", self$outcome_node)
      return(param_form)
    },
    cf_likelihood = function() {
      return(self$observed_likelihood)
    },
    update_nodes = function() {
      return(NULL)
    }
  ),
  private = list(
    .type = "E(Y)"
  )
)
