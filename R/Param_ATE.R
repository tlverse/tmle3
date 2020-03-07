#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
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
Param_ATE <- R6Class(
  classname = "Param_ATE",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node)
      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list_treatment <- c(intervention_list_treatment, censoring_intervention)
        intervention_list_control <- c(intervention_list_control, censoring_intervention)
      }

      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      HA_treatment <- cf_pA_treatment / pA
      HA_control <- cf_pA_control / pA

      # collapse across multiple intervention nodes
      if (!is.null(ncol(HA_treatment)) && ncol(HA_treatment) > 1) {
        HA_treatment <- apply(HA_treatment, 1, prod)
      }

      # collapse across multiple intervention nodes
      if (!is.null(ncol(HA_control)) && ncol(HA_control) > 1) {
        HA_control <- apply(HA_control, 1, prod)
      }

      HA <- HA_treatment - HA_control

      HA <- bound(HA, c(-40, 40))
      return(list(Y = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]


      # todo: make sure we support updating these params
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      Y <- tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE)

      EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
      EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)

      psi <- mean(EY1 - EY0)

      IC <- HA * (Y - EY) + (EY1 - EY0) - psi

      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
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
      return(c(self$outcome_node))
    }
  ),
  private = list(
    .type = "ATE",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE
  )
)
