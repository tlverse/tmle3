#' Define Likelihoods for Static Interventions
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_static <- R6Class(
  classname = "LF_static",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, value, ...) {
      private$.name <- name
      private$.value <- value
      private$.variable_type <- variable_type("constant", value)
    },
    get_likelihood = function(task, only_observed = FALSE) {
      
        observed <- task$get_tmle_node(self$name)
        likelihood <- as.numeric(self$value == observed)
      
      return(likelihood)
    }
  ),
  active = list(
    value = function() {
      return(private$.value)
    }
  ),
  private = list(
    .name = NULL,
    .value = NULL,
    .is_degenerate = TRUE
  )
)
