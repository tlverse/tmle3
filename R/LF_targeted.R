#' Use a likelihood factor from an existing targeted likelihood
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_fit, name, learner, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{learner}}{An sl3 learner to be used to estimate the factor
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{learner}}{The learner or learner fit object}
#'     }
#'
#' @export
LF_targeted <- R6Class(
  classname = "LF_targeted",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, base_likelihood, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.base_likelihood <- base_likelihood
    },
    get_likelihood = function(tmle_task, fold_number) {
      lik <- self$base_likelihood$get_likelihood(tmle_task, self$name, fold_number)

      return(lik)
    }
  ),
  active = list(
    base_likelihood = function() {
      return(private$.base_likelihood)
    }
  ),
  private = list(
    .name = NULL,
    .base_likelihood = NULL
  )
)
