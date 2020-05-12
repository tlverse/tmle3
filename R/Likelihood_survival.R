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
Likelihood_survival <- R6Class(
  classname = "Likelihood_survival",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  public = list(
    initialize = function(factor_list, cache = NULL, ...) {
      super$initialize(factor_list, cache, ...)
    },
    # TODO: get likelihoods for W, A
    get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full") {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      # TODO: decide which nodes to get
      nodes <- nodes[1:2]
      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number)
        })
        likelihood_dt <- as.data.table(all_likelihoods)
        setnames(likelihood_dt, nodes)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
      }
    },
    get_hazards = function(tmle_task, nodes = NULL, fold_number = "full") {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      # TODO: decide which nodes to get
      nodes <- nodes[3:4]
      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number)
        })
        likelihood_dt <- as.data.table(all_likelihoods)

        # TODO: check
        new_node_names <- unlist(lapply(nodes, convert_node_name, 
          ncol(likelihood_dt) / length(nodes)))
        setnames(likelihood_dt, new_node_names)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
      }
    },
    get_survival = function(tmle_task, nodes = NULL, fold_number = "full") {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      # TODO: decide which nodes to get
      nodes <- nodes[3:4]
      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number)
        })
        # TODO: transform hazards into survival
        all_likelihoods_surv <- lapply(seq(length(all_likelihoods)), function(i) {
          t(apply(1 - all_likelihoods[[i]], 1, cumprod))
        })
        likelihood_dt <- as.data.table(all_likelihoods_surv)

        # TODO: rename
        survival_names <- c("S_N", "S_A_c")
        new_node_names <- unlist(lapply(survival_names, convert_node_name, 
          ncol(likelihood_dt) / length(survival_names)))
        setnames(likelihood_dt, new_node_names)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
      }
    }
  ),
  active = list(),
  private = list()
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname Likelihood_survival
#'
#' @export
#
make_Likelihood_survival <- Likelihood_survival$new
