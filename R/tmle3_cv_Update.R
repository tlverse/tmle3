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
      # so we estimate epsilon using validation sets

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
    check_convergence = function(tmle_task) {
      ED_criterion <- 1 / tmle_task$nrow
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, cv_fold = 0)
        }
      )
      ICs <- sapply(estimates, `[[`, "IC")
      ED <- colMeans(ICs)
      return(max(abs(ED)) < ED_criterion)
    }
  )
)
