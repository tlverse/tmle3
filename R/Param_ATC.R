#' Additive Effect of Treatment Among the Treated
#'
#' Parameter definition for the Additive Effect of Treatment Among the Treated (ATT). Currently supports multiple static intervention nodes.
#' Does yet not support dynamic rule or stochastic interventions.
#'
#' @section Current Issues:
#' \itemize{
#'   \item clever covariates doesn't support updates; always uses initial (necessary for iterative TMLE, e.g. stochastic intervention)
#'   \item doesn't integrate over possible counterfactuals (necessary for stochastic intervention)
#'   \item clever covariate gets recalculated all the time (inefficient)
#' }
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
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_ATC <- R6Class(
  classname = "Param_ATC",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {

      # flip treatment and control
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      private$.outcome_node <- outcome_node
      private$.param_att <- Param_ATT$new(observed_likelihood, intervention_list_control, intervention_list_treatment, outcome_node)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      att_cc <- self$param_att$clever_covariates(tmle_task, fold_number)

      atc_cc <- list(A = -1 * att_cc$A, Y = -1 * att_cc$Y)
      return(atc_cc)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      att_est <- self$param_att$estimates(tmle_task, fold_number)
      result <- list(psi = -1 * att_est$psi, IC = -1 * att_est$IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATC[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    update_nodes = function() {
      return(c(self$outcome_node, names(self$intervention_list_treatment)))
    },
    param_att = function() {
      return(private$.param_att)
    }
  ),
  private = list(
    .type = "ATC",
    .param_att = NULL,
    .outcome_node = NULL,
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL
  )
)
