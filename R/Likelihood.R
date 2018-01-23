#' Class for Likelihood Computations
#'
# TODO: different factor types
#   marginal - empirical density of W Q
#   conditional pdf - intervention nodes g's
#   conditional mean - outcome node Q_bar
#   integrate conditional mean Q_bar wrt to intervention g's.
#   If degenerate, this means just plugging in intervention values, otherwise
#   have to actually do discrete expectation
#
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
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
    get_factor = function(factor_name) {
      return(self$factor_list[[factor_name]])
    },
    modify_factors = function(new_factor_list) {
      new_factor_names <- sapply(new_factor_list, `[[`, "name")
      factor_list <- self$factor_list
      factor_list[new_factor_names] <- new_factor_list

      new_likelihood <- self$clone()
      new_likelihood$initialize(factor_list)
      new_likelihood$set_train(self$fit_object, self$training_task)
      return(new_likelihood)
    },
    validate_task = function(task) {
      assert_that(is(task, "tmle_Task"))

      factor_list <- self$factor_list
      factor_names <- names(factor_list)
      task_nodes <- names(task$tmle_nodes)
      if (!setequal(task_nodes, factor_names)) {
        stop("factor_list and task$tmle_nodes must have matching names")
      }
    },
    get_likelihoods = function(task, nodes = NULL) {
      self$validate_task(task)
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }
      likelist <- list()
      for (likelihood_factor in factor_list) {
        factor_name <- likelihood_factor$name
        likes <- likelihood_factor$get_likelihood(task, only_observed = TRUE)
        likelist[[factor_name]] <- likes
      }
      if (length(likelist) > 1) {
        return(as.data.table(likelist))
      } else {
        return(likelist[[1]])
      }
    },
    joint_likelihoods = function(task, nodes=NULL){
      likelihoods <- self$get_likelihoods(task, nodes=nodes)
      
      #todo: check this works if nodes is length 1
      joint <- apply(likelihoods, 1, prod)
      
      return(joint)
    },
    get_possible_counterfacutals = function(task, nodes=NULL){
      
      #get factors for nodes
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }
      
      all_levels <- lapply(factor_list, function(likelihood_factor){likelihood_factor$variable_type$levels})
      all_levels <- all_levels[ !(sapply(all_levels, is.null))]
      level_grid <- expand.grid(all_levels)
      return(level_grid)
    },
    E_f_x = function(tmle_task, f_x){
      cf_grid <- self$get_possible_counterfacutals(task)
      # todo: rewrite this so it does't recalculate likelihoods (ie recursively)
      # should start with base node and go forward
      prods <- lapply(seq_len(nrow(cf_grid)), function(cf_row){
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), as.data.table(cf_grid[cf_row,]))
        likelihoods <- self$joint_likelihoods(cf_task)
        f_vals <- f_x(cf_task)
        return(f_vals*likelihoods)
      })
      
      prodmat = do.call(cbind, prods)
      result <- sum(prodmat)
      
      return(result)
    },
    
    EY = function(tmle_task, mean_node_name="Y"){
      # identify set of all ancestors
      nodes <- tmle_task$tmle_nodes
      mean_node <- nodes[[mean_node_name]]
      ancestor_nodes <- all_ancestors(mean_node_name, nodes)
      mean_factor <- self$factor_list[[mean_node_name]]
      # get cf possibilities only for these ancestors
      cf_grid <- self$get_possible_counterfacutals(tmle_task, ancestor_nodes)
      
      prods <- lapply(seq_len(nrow(cf_grid)), function(cf_row){
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), as.data.table(cf_grid[cf_row,, drop=FALSE]))
        ey <- mean_factor$get_prediction(cf_task)
        likelihoods <- self$joint_likelihoods(cf_task, ancestor_nodes)
        return(ey*likelihoods)
      })
      
      prodmat = do.call(cbind, prods)
      result <- sum(prodmat)
      
      return(result)
    },
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
      fit_object <- private$.train(task)
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
    }
  ),
  private = list(
    .train = function(tmle_task) {
      # TODO: move some of this to .pretrain so we can delay it
      for (likelihood_factor in self$factor_list) {
        likelihood_factor$train(tmle_task)
      }
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do
      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood")
    }
  )
)
