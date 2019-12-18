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
Param_ATT <- R6Class(
  classname = "Param_ATT",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # todo: actually the union of the treatment and control nodes?
      intervention_nodes <- names(self$intervention_list_treatment)

      # todo: make sure we support updating these params
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      cf_task_treatment <- self$cf_likelihood_treatment$cf_tasks[[1]]
      cf_task_control <- self$cf_likelihood_control$cf_tasks[[1]]

      pA1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, intervention_nodes, fold_number)
      pA1_overall <- mean(pA1)

      HA <- (cf_pA_treatment - cf_pA_control * (pA1 / (1 - pA1)))



      EY1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, self$outcome_node, fold_number)
      EY0 <- self$observed_likelihood$get_likelihoods(cf_task_control, self$outcome_node, fold_number)

      psi <- mean((EY1 - EY0) * (pA1 / pA1_overall))
      CY <- (EY1 - EY0) - psi

      return(list(A = CY, Y = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # todo: actually the union of the treatment and control nodes?
      intervention_nodes <- names(self$intervention_list_treatment)

      # todo: make sure we support updating these params
      # pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # pA_overall <- mean(pA)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)



      cf_task_treatment <- self$cf_likelihood_treatment$cf_tasks[[1]]
      cf_task_control <- self$cf_likelihood_control$cf_tasks[[1]]

      pA1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, intervention_nodes, fold_number)
      pA1_overall <- mean(pA1)

      EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
      EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)

      psi <- mean((EY1 - EY0) * (pA1 / pA1_overall))

      Y <- tmle_task$get_tmle_node(self$outcome_node)

      clever_covariates <- self$clever_covariates(tmle_task, fold_number)
      HA <- clever_covariates$Y
      CY <- clever_covariates$A

      IC <- (HA * (Y - EY) + CY * cf_pA_treatment) / pA1_overall

      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATT[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
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
    }
  ),
  private = list(
    .type = "ATT",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL
  )
)
