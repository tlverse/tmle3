#' Counterfactual likelihood
#'
#' Based on a sample_likelihood estimated from observed data.
#' Modifies likelihood factors according to intervention_list
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#'
#' @export
CF_likelihood <- R6Class(
  classname = "CF_likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  public = list(
    initialize = function(sample_likelihood, intervention_list, ...) {
      private$.sample_likelihood <- sample_likelihood
      private$.intervention_list <- intervention_list
      params <- args_to_list()
      super$initialize(params)
    }
  ),
  active = list(
    name = function() {
      node_names <- names(self$intervention_list)
      node_values <- sapply(self$intervention_list, `[[`, "values")
      intervention_name <- paste(sprintf("%s=%s", node_names, as.character(node_values)), collapse = ", ")
      return(intervention_name)
    },
    sample_likelihood = function() {
      return(private$.sample_likelihood)
    },
    intervention_list = function() {
      return(private$.intervention_list)
    },
    factor_list = function() {
      fl <- self$sample_likelihood$factor_list
      fl[names(self$intervention_list)] <- self$intervention_list[names(self$intervention_list)]
      return(fl)
    },
    update_list = function() {
      self$sample_likelihood$update_list
    }
  ),
  private = list(
    .sample_likelihood = NULL,
    .intervention_list = NULL
  )
)
