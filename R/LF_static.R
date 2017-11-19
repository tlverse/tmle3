#' Defining Likelihood Functionals for Static Interventions
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_static <- R6Class(classname = "LF_static",
                     portable = TRUE,
                     class = TRUE,
                     inherit = LF_base,
  public = list(
    initialize <- function(name, value, ...) {
      private$.name <- name
      private$.value <- value
    },
    get_values <- function(task) {
      values <- rep(self$value, task$nrow)
      return(values)
    },
    get_likelihood <- function(task, only_observed = FALSE) {
      node_task <- task$get_regression_task(self$name)
      values <- self$get_values(task)
      outcome_type <- node_task$outcome_type

      if (only_observed) {
        observed <- outcome_type$format(node_task$Y)
        likelihood <- as.numeric(values == observed)
      } else {
        if (outcome_type$type == "binomial") {
          levels <- outcome_type$levels
          level_mat <- matrix(levels, nrow = length(values),
                              ncol = length(levels), byrow = TRUE)
          likelihood <- apply(level_mat, MARGIN = 2,
                              function(level_vec) {
                                as.numeric(level_vec == values)
                              })
        } else {
          stop("currently, only binomial likelihoods are supported")
        }
      }
      return(likelihood)
    }
  ),
  active = list(
    value <- function() {
      return(private$.value)
    }
  ),
  private = list(
    .name = NULL,
    .value = NULL,
    .is_degenerate = TRUE
  )
)

