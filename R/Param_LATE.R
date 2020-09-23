#' Longitudinal Average Treatment Effect
#' Parameter for ATE in (L0 A0 L1 A1 ... Y) data structure with any number of time dependent covariates and treatments.
#' Supports arbitrarily many time points.
#' Parameter definition for the Longitudinal Average Treatment Effect (LATE).
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
Param_LATE <- R6Class(
  classname = "Param_LATE",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      all_nodes <- names(observed_likelihood$training_task$npsem)
      A_nodes <- grep("A", all_nodes, value = T)
      L_nodes <- grep("L", all_nodes, value = T)
      private$.update_nodes <- c(L_nodes, outcome_node)

      super$initialize(observed_likelihood, list(), outcome_node)



      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      # Train the gradient
      private$.gradient <- Gradient$new(observed_likelihood,
                                        ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment, cf_likelihood_control = self$cf_likelihood_control),
                                        projection_task_generator = gradient_generator_late,
                                        target_nodes =  self$update_nodes)

      if(inherits(observed_likelihood, "Targeted_Likelihood")){
        fold_number <- observed_likelihood$updater$update_fold
      } else {
        fold_number <- "full"
      }
      private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      print(tmle_task$uuid)
      print(node)
      update_nodes <- intersect(self$update_nodes, attr(tmle_task, "target_nodes"))
      if(!is.null(node)){
        update_nodes <- c(node)
      }
      islong = F
      if(is.null(update_nodes)){
        update_nodes <- self$update_nodes
      } else {
      islong= T
      }
      print(update_nodes)
      EICs <- lapply(update_nodes, function(node){
        return(self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC)
      })

      names(EICs) <- update_nodes
      return(EICs)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIC <- (do.call(cbind, self$clever_covariates(tmle_task, fold_number)))

      #TODO need to montecarlo simulate from likleihood to eval parameter.

      # todo: make sure we support updating these params
      # pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      #
      # # todo: extend for stochastic
      # cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      # cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      #
      # Y <- tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE)
      #
      # EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      # EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
      # EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)
      #
      # psi <- mean(EY1 - EY0)
      #
      # IC <- EIC + (EY1 - EY0) - psi

      psi = rep(0, length(EIC))
      IC <- rowSums(EIC)
      result <- list(psi = psi, IC = IC, EIC = colMeans(EIC))
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
      return(c(private$.update_nodes))
    },
    gradient = function(){
      private$.gradient
    }
  ),
  private = list(
    .type = "ATE",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = FALSE,
    .gradient = NULL,
    .submodel_type_supported = c("EIC"),
    .update_nodes = NULL
  )
)

