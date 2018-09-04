#' Defines an update (submodel+loss function)
#'
#' Current Limitations:
#' only updating one node
#' loss function and submodel are hard-coded (need to accept arguments for these)
#' no support for one-step (recursive TMLE)
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_cv_Update <- R6Class(
  classname = "tmle3_cv_Update",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Update,
  public = list(
    update_step = function(likelihood, tmle_task, cv_fold = 0) {
      # cv_fold=0 -- validation sets
      # so we estimate epsilon using valudation sets

      # get new submodel fit
      all_submodels <- self$generate_submodel_data(likelihood, tmle_task, cv_fold)
      new_epsilons <- self$fit_submodels(all_submodels)

      # update likelihoods
      # todo: think more carefully about what folds to update
      likelihood$update(new_epsilons, self$step_number, -1)
      likelihood$update(new_epsilons, self$step_number, 0)

      # increment step count
      private$.step_number <- private$.step_number + 1
    },
    generate_submodel_data = function(likelihood, tmle_task, cv_fold = -1) {
      update_nodes <- self$update_nodes

      # todo: support not getting observed for case where we're applying updates instead of fitting them
      clever_covariates <- lapply(self$tmle_params, function(tmle_param) tmle_param$clever_covariates(tmle_task, cv_fold))

      observed_values <- lapply(update_nodes, tmle_task$get_tmle_node, bound = TRUE)

      all_submodels <- lapply(update_nodes, function(update_node) {
        node_covariates <- lapply(clever_covariates, `[[`, update_node)
        covariates_dt <- do.call(cbind, node_covariates)
        observed <- tmle_task$get_tmle_node(update_node, bound = TRUE)
        initial <- likelihood$get_likelihood(tmle_task, update_node, cv_fold)
        submodel_data <- list(
          observed = observed,
          H = covariates_dt,
          initial = initial
        )
      })

      names(all_submodels) <- update_nodes

      return(all_submodels)
    },
    check_convergence = function(tmle_task){
      ED_criterion <- 1 / tmle_task$nrow
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, cv_fold=0)
        }
      )
      ICs <- sapply(estimates, `[[`, "IC")
      ED <- colMeans(ICs)
      return(max(abs(ED)) < ED_criterion)
    }
  )
)
