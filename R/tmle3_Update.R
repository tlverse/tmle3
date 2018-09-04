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
tmle3_Update <- R6Class(
  classname = "tmle3_Update",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(maxit=100) {
      private$.maxit=maxit
    },
    update_step = function(likelihood, tmle_task) {


      # get new submodel fit
      all_submodels <- self$generate_submodel_data(likelihood, tmle_task)
      new_epsilons <- self$fit_submodels(all_submodels)

      # update likelihoods
      likelihood$update(new_epsilons, self$step_number)

      # increment step count
      private$.step_number <- private$.step_number + 1
    },
    generate_submodel_data = function(likelihood, tmle_task) {
      update_nodes <- self$update_nodes

      # todo: support not getting observed for case where we're applying updates instead of fitting them
      clever_covariates <- lapply(self$tmle_params, function(tmle_param) tmle_param$clever_covariates(tmle_task))

      observed_values <- lapply(update_nodes, tmle_task$get_tmle_node, bound = TRUE)

      all_submodels <- lapply(update_nodes, function(update_node) {
        node_covariates <- lapply(clever_covariates, `[[`, update_node)
        covariates_dt <- do.call(cbind, node_covariates)
        observed <- tmle_task$get_tmle_node(update_node, bound = TRUE)
        initial <- likelihood$get_likelihood(tmle_task, update_node)
        submodel_data <- list(
          observed = observed,
          H = covariates_dt,
          initial = initial
        )
      })

      names(all_submodels) <- update_nodes

      return(all_submodels)
    },
    fit_submodel = function(submodel_data) {
      suppressWarnings({
        submodel_fit <- glm(observed ~ H - 1, submodel_data, offset = qlogis(submodel_data$initial), family = binomial())
      })
      epsilon <- coef(submodel_fit)

      # this protects against collinear covariates (which we don't care about, we just want an update)
      epsilon[is.na(epsilon)] <- 0
      return(epsilon)
    },
    fit_submodels = function(all_submodels) {
      all_epsilon <- lapply(all_submodels, self$fit_submodel)

      names(all_epsilon) <- names(all_submodels)
      private$.epsilons <- c(private$.epsilons, list(all_epsilon))

      return(all_epsilon)
    },
    submodel = function(epsilon, initial, H) {
      plogis(qlogis(initial) + H %*% epsilon)
    },
    loss_function = function(estimate, observed) {
      -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
    },
    apply_submodel = function(submodel_data, epsilon) {
      self$submodel(epsilon, submodel_data$initial, submodel_data$H)
    },
    apply_submodels = function(all_submodels, all_epsilon) {
      update_nodes <- self$update_nodes
      updated_likelihood <- mapply(self$apply_submodel, all_submodels, all_epsilon, SIMPLIFY = FALSE)
      return(updated_likelihood)
    },
    check_convergence = function(tmle_task){
      ED_criterion <- 1 / tmle_task$nrow
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task)
        }
      )
      ICs <- sapply(estimates, `[[`, "IC")
      ED <- colMeans(ICs)
      return(max(abs(ED)) < ED_criterion)
    },
    update = function(likelihood, tmle_task){
      
      maxit <- private$.maxit
      for (steps in seq_len(maxit)) {
        self$update_step(likelihood, tmle_task)
        if(self$check_convergence(tmle_task)){
          break
        }
      }
      
    }
  ),
  active = list(
    epsilons = function() {
      return(private$.epsilons)
    },
    tmle_params = function(new_params = NULL) {
      if (!is.null(new_params)) {
        if (inherits(new_params, "Param_base")) {
          new_params <- list(new_params)
        }
        private$.tmle_params <- new_params
        private$.update_nodes <- unique(unlist(lapply(new_params, `[[`, "update_nodes")))
      }
      return(private$.tmle_params)
    },
    update_nodes = function() {
      return(private$.update_nodes)
    },
    step_number = function() {
      return(private$.step_number)
    }
  ),
  private = list(
    .epsilons = list(),
    .tmle_params = NULL,
    .update_nodes = NULL,
    .step_number = 0,
    .maxit = 100
  )
)
