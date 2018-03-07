#' Likelihood Factor Estimated using Empirical Distribution
#'
#' Uses the empirical probability distribution (puts mass \eqn{1/n} on each of the observations, or uses weights if specified) to estimate a marginal density.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#' Only compatible with marginal likelihoods (no parent nodes). Only compatible with densities (no conditional means).
#' The \code{type} argument will be ignored if specified.
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
#'   \code{define_lf(LF_emp, name, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @export
LF_emp <- R6Class(
  classname = "Lf_np",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, ...) {
      super$initialize(name, ..., type = "density")
      private$.name <- name
    },
    get_mean = function(tmle_task) {
      stop("nothing to predict")
    },
    get_likelihood = function(tmle_task) {
      weights <- tmle_task$weights
      return(weights / sum(weights))
    }
  ),
  active = list(),
  private = list(
    .name = NULL
  )
)
