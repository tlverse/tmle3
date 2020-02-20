#' Stratified Parameter Estimates
#'
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
#'   \code{define_param(Param_TSM, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
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
Param_stratified <- R6Class(
  classname = "Param_stratified",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, param_base, strata_variable, ..., outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.param_base <- param_base
      private$.type <- sprintf("stratified %s", param_base$type)
      private$.strata_variable <- strata_variable

      V <- observed_likelihood$training_task$get_data(, strata_variable)

      strata <- V[, list(weight = observed_likelihood$training_task$nrow / .N), by = names(V)]
      set(strata, , "strata_i", factor(1:nrow(strata)))
      private$.strata <- strata
    },
    
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      base_covs <- self$param_base$clever_covariates(tmle_task, fold_number)
      strata_weights <- self$get_strata_weights(tmle_task)

      strata_covs <- lapply(base_covs, `*`, strata_weights)
      return(strata_covs)
    },
    
    estimates = function(tmle_task = NULL, fold_number = "full") {
      strata_weights <- self$get_strata_weights(tmle_task)
      strata_tasks <- apply(strata_weights, 2, function(weights) tmle_task[which(weights != 0)])
      strata_ests <- lapply(strata_tasks, self$param_base$estimates, fold_number)
      psi <- sapply(strata_ests, `[[`, "psi")

      IC <- strata_weights
      all_ICs <- unlist(lapply(strata_ests, `[[`, "IC"))
      
      IC[which(strata_weights != 0)] <- IC[which(strata_weights != 0)] * all_ICs
      result <- list(psi = psi, IC = IC)
      return(result)
    },
    
    get_strata_weights = function(tmle_task) {
      V <- tmle_task$get_data(, self$strata_variable)
      strata <- self$strata
      combined <- merge(V, strata, by = self$strata_variable, sort = FALSE, all.x = TRUE)
      combined[, index := .I]
      strata_weights_dt <- dcast(combined, index ~ strata_i, value.var = "weight", fill = 0, drop = FALSE)
      strata_weights <- as.matrix(strata_weights_dt[, -1, with = FALSE])
      return(strata_weights)
    },
    
    get_strata_indicators = function(tmle_task) {
      # deprecated
      # keeping it here temporarily in case it's needed
      V <- tmle_task$get_data(, self$strata_variable)
      strata <- self$strata[, !"weight"][, "indicators":= 1]
      combined <- merge(V, strata, by = self$strata_variable, sort = FALSE, all.x = TRUE)
      combined[, index := .I]
      strata_indicators_dt <- dcast(combined, index ~ strata_i, value.var = "indicators", fill = 0, drop = FALSE)
      strata_indicators <- as.matrix(strata_indicators_dt[, -1, with = FALSE])
      return(strata_indicators)
    }
    
  ),
  active = list(
    name = function() {
      strata_labels <- apply(self$strata[, self$strata_variable, with = FALSE], 1, paste, collapse = ", ")
      param_form <- sprintf(
        "%s | V=%s",
        self$param_base$name,
        strata_labels
      )
      return(param_form)
    },
    param_base = function() {
      return(private$.param_base)
    },
    strata_variable = function() {
      return(private$.strata_variable)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    strata = function() {
      return(private$.strata)
    }
  ),
  private = list(
    .type = NULL,
    .param_base = NULL,
    .strata_variable = NULL,
    .strata = NULL,
    .supports_outcome_censoring = TRUE
  )
)
