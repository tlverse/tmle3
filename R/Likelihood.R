#' Class for Likelihood Computations
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom delayed delayed_fun bundle_delayed
#'
#' @export
#
Likelihood <- R6Class(
  classname = "Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(factor_list, ...) {
      params <- args_to_list()

      factor_names <- sapply(factor_list, `[[`, "name")
      names(factor_list) <- factor_names
      params$factor_list <- factor_list

      super$initialize(params)
    },
    print = function() {
      lapply(self$factor_list, print)
      invisible(NULL)
    },
    validate_task = function(tmle_task) {
      assert_that(is(tmle_task, "tmle3_Task"))

      factor_list <- self$factor_list
      factor_names <- names(factor_list)
      task_nodes <- names(tmle_task$tmle_nodes)
      if (!setequal(task_nodes, factor_names)) {
        stop("factor_list and task$tmle_nodes must have matching names")
      }
    },
    get_initial_likelihoods = function(tmle_task, nodes = NULL, only_observed = TRUE) {
      self$validate_task(tmle_task)
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }
      likelist <- list()
      for (likelihood_factor in factor_list) {
        factor_name <- likelihood_factor$name
        if (likelihood_factor$type == "mean") {
          likes <- likelihood_factor$get_prediction(tmle_task)
        } else {
          likes <- likelihood_factor$get_likelihood(tmle_task, only_observed = only_observed)
        }
        likelist[[factor_name]] <- likes
      }
      # if (length(likelist) > 1) {
      return(as.data.table(likelist))
      # } else {
      # return(likelist[[1]])
      # }
    },
    get_likelihoods = function(tmle_task, nodes = NULL, only_observed = TRUE) {
      # todo: maybe get all likelihoods here, then subset before returning
      initial_likelihoods <- self$get_initial_likelihoods(tmle_task, nodes = nodes, only_observed = only_observed)

      # apply updates
      if (is.null(self$update_list)) {
        return(initial_likelihoods)
      } else {
        updater <- self$update_list
        updater$epsilons
        updated_likeilhoods <- updater$apply_updates(tmle_task, self, initial_likelihoods)
        return(updated_likeilhoods)
      }
    },
    joint_likelihoods = function(task, nodes=NULL) {
      likelihoods <- self$get_likelihoods(task, nodes = nodes)

      # todo: check this works if nodes is length 1
      joint <- apply(likelihoods, 1, prod)

      return(joint)
    },
    get_possible_counterfacutals = function(nodes=NULL) {

      # get factors for nodes
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }

      all_levels <- lapply(factor_list, function(likelihood_factor) {
        likelihood_factor$variable_type$levels
      })
      all_levels <- all_levels[ !(sapply(all_levels, is.null))]
      level_grid <- expand.grid(all_levels)
      return(level_grid)
    },
    # E_f_x = function(tmle_task, f_x) {
    #   cf_grid <- self$get_possible_counterfacutals()
    #   # todo: rewrite this so it does't recalculate likelihoods (ie recursively)
    #   # should start with base node and go forward
    #   prods <- lapply(seq_len(nrow(cf_grid)), function(cf_row) {
    #     cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), as.data.table(cf_grid[cf_row, ]))
    #     likelihoods <- self$joint_likelihoods(cf_task)
    #     f_vals <- f_x(cf_task)
    #     return(f_vals * likelihoods)
    #   })
    # 
    #   prodmat <- do.call(cbind, prods)
    #   result <- sum(prodmat)
    # 
    #   return(result)
    # },
    # EY = function(tmle_task, mean_node_name="Y") {
    #   # identify set of all ancestors
    #   nodes <- tmle_task$tmle_nodes
    #   mean_node <- nodes[[mean_node_name]]
    #   ancestor_nodes <- all_ancestors(mean_node_name, nodes)
    #   mean_factor <- self$factor_list[[mean_node_name]]
    #   # get cf possibilities only for these ancestors
    #   cf_grid <- self$get_possible_counterfacutals(ancestor_nodes)
    # 
    #   prods <- lapply(seq_len(nrow(cf_grid)), function(cf_row) {
    #     cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), as.data.table(cf_grid[cf_row, , drop = FALSE]))
    #     ey <- mean_factor$get_prediction(cf_task)
    #     likelihoods <- self$joint_likelihoods(cf_task, ancestor_nodes)
    #     return(ey * likelihoods)
    #   })
    # 
    #   prodmat <- do.call(cbind, prods)
    #   result <- sum(prodmat)
    # 
    #   return(result)
    # },
    get_predictions = function(task, nodes = NULL) {
      self$validate_task(task)
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }
      predlist <- list()
      for (likelihood_factor in factor_list) {
        if (inherits(likelihood_factor, "LF_fit")) {
          factor_name <- likelihood_factor$name
          preds <- likelihood_factor$get_prediction(task)
          predlist[[factor_name]] <- preds
        }
      }
      if (length(predlist) > 1) {
        return(as.data.table(predlist))
      } else {
        return(predlist[[1]])
      }
    },
    base_train = function(task, pretrain) {
      self$validate_task(task)
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    },
    base_predict = function(task = NULL) {
      if (is.null(task)) {
        task <- private$.training_task
      }
      self$validate_task(task)
      predictions <- private$.predict(task)
      return(predictions)
    },
    base_chain = function(task = NULL) {
      if (is.null(task)) {
        task <- private$.training_task
      }
      self$validate_task(task)
      predictions <- private$.chain(task)
      return(predictions)
    }
  ),
  active = list(
    factor_list = function() {
      return(self$params$factor_list)
    },
    observed_values = function() {
      if (is.null(private$.observed_values)) {
        private$.observed_values <- self$get_likelihoods(self$training_task)
      }
      return(private$.observed_values)
    },
    update_list = function(new_update_list = NULL) {
      if (!is.null(new_update_list)) {
        private$.update_list <- new_update_list
      }
      return(private$.update_list)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task) {
      # TODO: move some of this to .pretrain so we can delay it
      #don't nest delayed calls
      #get the delayed calls for training
      #then pass them to other delayed methods
      #remember that tmle_task could be delayed too
      #delay get regression task
      #get learner
      #delay_learner_train the results
      #set in .train via setter
      #how do to this but still respect train methods for LFs?

      factor_fits <- lapply(self$factor_list, function(factor) factor$delayed_train(tmle_task))
      result <- bundle_delayed(factor_fits)
      return(result)

    },
    .train = function(tmle_task, factor_fits) {
      factor_list <- self$factor_list
      for (i in seq_along(factor_list)) {
          factor_list[[i]]$train(tmle_task, factor_fits[[i]])
      }
      #set learners here
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do
      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood. See Likelihood$get_likelihoods for analogous method")
    },
    .observed_values = NULL,
    .update_list = NULL
  )
)
