#' Class for Estimation of TSM Parameter
#'
#' Current Limitations:
#' clever covariates doesn't support updates; always uses initial (necessary for stochastic intervention)
#' doesn't integrate over possible counterfactuals (necessary for stochastic intervention)
#' clever covariate gets recalculated all the time (inefficient)
#' @importFrom R6 R6Class
#'
#' @export
#
Param_TSM <- R6Class(
  classname = "Param_TSM",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(intervention_list, outcome_node = "Y") {
      private$.intervention_list <- intervention_list
      private$.outcome_node <- outcome_node
    },
    setup = function(sample_likelihood) {
      # todo: this is something like a fit operation
      private$.sample_likelihood <- sample_likelihood
      private$.cf_likelihood <- CF_likelihood$new(sample_likelihood, self$intervention_list)
    },
    clever_covariates = function(tmle_task) {
      intervention_nodes <- names(self$intervention_list)
      # todo: make sure we support updating these params
      pA <- self$sample_likelihood$get_initial_likelihoods(tmle_task, intervention_nodes)
      cf_pA <- self$cf_likelihood$get_initial_likelihoods(tmle_task, intervention_nodes)

      HA <- cf_pA / pA

      # collapse across multiple intervention nodes
      if (!is.null(ncol(HA)) && ncol(HA) > 1) {
        HA <- apply(HA, 1, prod)
      }
      return(list(Y = unlist(HA, use.names=FALSE)))
    },
    estimates = function(tmle_task) {
      # todo: extend for stochastic interventions
      cfs <- self$cf_likelihood$get_possible_counterfacutals()
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cfs)

      
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      
      
      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task)[[self$outcome_node]]

      # clever_covariates happen here (for all fit params), but this is repeated computation
      EY <- unlist(self$sample_likelihood$get_likelihoods(tmle_task, self$outcome_node), use.names=FALSE)

      # clever_covariates happen here (for all fit params)
      EY1 <- unlist(self$cf_likelihood$get_likelihoods(cf_task, self$outcome_node), use.names=FALSE)
      
      # todo: integrate unbounding logic into likelihood class, or at least put it in a function
      variable_type <- tmle_task$tmle_nodes[[self$outcome_node]]$variable_type
      if((variable_type$type=="continuous")&&(!is.na(variable_type$bounds))){
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
    intervention_list = function() {
      return(private$.intervention_list)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    sample_likelihood = function() {
      return(private$.sample_likelihood)
    },
    update_nodes = function() {
      return(self$outcome_node)
    }
  ),
  private = list(
    .intervention_list = NULL,
    .sample_likelihood = NULL,
    .cf_likelihood = NULL
  )
)
