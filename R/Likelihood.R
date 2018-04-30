#' Class for Likelihood
#'
#' This object represents an estimate of the relevant factors of the likelihood estimated from data, or based on \emph{a priori} knowledge where appropriate.
#' That is, it represents some subset of $P_n$. This object inherits from \code{\link[sl3]{Lrnr_base}}, and so shares some properties with \code{sl3} learners.
#' Specifically, to fit a likelihood object to data, one calls \code{likelihood$train(tmle3_task)}.
#' Each likelihood factor is represented by an object inheriting from \code{\link{LF_base}}.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base
#' @importFrom assertthat assert_that is.count is.flag
#' @importFrom delayed bundle_delayed
#' @import data.table
#' @family Likelihood objects
#' @export
#'
#' @keywords data
#'
#' @return \code{Likelihood} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Likelihood_extra
#'
#' @export
Likelihood <- R6Class(
  classname = "Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(factor_list, cache = NULL, ...) {
      params <- args_to_list()
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }

      factor_names <- sapply(factor_list, `[[`, "name")
      names(factor_list) <- factor_names
      params$factor_list <- factor_list
      if (is.null(cache)) {
        cache <- Likelihood_cache$new()
      }
      private$.cache <- cache

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
      task_nodes <- names(tmle_task$npsem)
      if (!all(factor_names %in% task_nodes)) {
        stop("factor_list and task$npsem must have matching names")
      }
    },
    get_likelihood = function(tmle_task, node) {
      likelihood_factor <- self$factor_list[[node]]
      # first check for cached values for this task
      likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task)

      if (is.null(likelihood_values)) {
        # if not, generate new ones
        likelihood_values <- likelihood_factor$get_likelihood(tmle_task)
        self$cache$set_values(likelihood_factor, tmle_task, 0, likelihood_values)
      }

      return(likelihood_values)
    },
    get_likelihoods = function(tmle_task, nodes=NULL) {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node)
        })
        likelihood_dt <- as.data.table(all_likelihoods)
        setnames(likelihood_dt, nodes)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]]))
      }
    },
    get_possible_counterfactuals = function(nodes=NULL) {

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
    base_train = function(task, pretrain) {
      self$validate_task(task)
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    }
  ),
  active = list(
    factor_list = function() {
      return(self$params$factor_list)
    },
    nodes = function() {
      return(names(self$factor_list))
    },
    cache = function() {
      return(private$.cache)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task) {
      factor_fits <- lapply(self$factor_list, function(factor) factor$delayed_train(tmle_task))
      result <- bundle_delayed(factor_fits)
      return(result)
    },
    .train = function(tmle_task, factor_fits) {
      factor_list <- self$factor_list
      for (i in seq_along(factor_list)) {
        factor_list[[i]]$train(tmle_task, factor_fits[[i]])
      }
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do

      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood. See Likelihood$get_likelihoods for analogous method")
    },
    .chain = function(tmle_task) {
      stop("chain method doesn't work for Likelihood. Currently, no analogous functionality")
    },
    .cache = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname Likelihood
#'
#' @export
#
make_Likelihood <- Likelihood$new
