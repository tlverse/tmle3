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
#'   \code{define_param(Param_ATT, observed_likelihood, outcome_node, intervention_list, ...)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
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
    initialize = function(observed_likelihood, outcome_node = "Y", intervention_list_treatment, intervention_list_control) {
      super$initialize(observed_likelihood, outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      # todo: actually the union of the treatment and control nodes?
      intervention_nodes <- names(self$intervention_list_treatment)
      
      # todo: make sure we support updating these params
      pA <- self$observed_likelihood$get_initial_likelihoods(tmle_task, intervention_nodes)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_initial_likelihoods(tmle_task, intervention_nodes)
      cf_pA_control <- self$cf_likelihood_control$get_initial_likelihoods(tmle_task, intervention_nodes)
      
      #todo: rethink that last term
      HA <- cf_pA_treatment - cf_pA_control * ((1-pA)/pA)
      
      #todo: rethink variable names
      
      # todo: extend for stochastic interventions
      cfs <- self$cf_likelihood_treatment$get_possible_counterfacutals()
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cfs)
      
      
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      
      
      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task)[[self$outcome_node]]
      
      # clever_covariates happen here (for all fit params), but this is repeated computation
      EY <- unlist(self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node), use.names = FALSE)
      
      # clever_covariates happen here (for all fit params)
      EY1 <- unlist(self$cf_likelihood$get_likelihoods(cf_task, self$outcome_node), use.names = FALSE)
      
      EY1 <- 
      # collapse across multiple intervention nodes
      if (!is.null(ncol(HA)) && ncol(HA) > 1) {
        HA <- apply(HA, 1, prod)
      }
      return(list(Y = unlist(HA, use.names = FALSE)))
    },
    estimates = function(tmle_task = NULL) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # todo: extend for stochastic interventions
      cfs <- self$cf_likelihood$get_possible_counterfacutals()
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cfs)


      Y <- tmle_task$get_tmle_node(self$outcome_node)


      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task)[[self$outcome_node]]

      # clever_covariates happen here (for all fit params), but this is repeated computation
      EY <- unlist(self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node), use.names = FALSE)

      # clever_covariates happen here (for all fit params)
      EY1 <- unlist(self$cf_likelihood$get_likelihoods(cf_task, self$outcome_node), use.names = FALSE)

      # todo: integrate unbounding logic into likelihood class, or at least put it in a function
      variable_type <- tmle_task$npsem[[self$outcome_node]]$variable_type
      if ((variable_type$type == "continuous") && (!is.na(variable_type$bounds))) {
        bounds <- variable_type$bounds
        scale <- bounds[2] - bounds[1]
        shift <- bounds[1]
        EY <- EY * scale + shift
        EY1 <- EY1 * scale + shift
      }

      # todo: separate out psi
      psi <- mean(EY1)
      IC <- HA * (Y - EY) + EY1 - psi

      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[%s_{%s}]", self$outcome_node, self$cf_likelihood$name)
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
      return(self$outcome_node)
    }
  ),
  private = list(
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL
  )
)
