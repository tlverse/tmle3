
stub_known <- function(task) {
  stop("mean or density function was not provided to LF_known and then values were requested")
}

#' Known True Likelihood Factor
#'
#' Incorporate existing knowledge about the likelihood
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
#'   \code{define_lf(LF_fit, name, mean_fun, density_fun, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{mean_fun}}{A function that takes a sl3 regression task and returns true conditional means
#'     }
#'     \item{\code{density_fun}}{A function that takes a sl3 regression task and returns true conditional densities
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @export
LF_known <- R6Class(
  classname = "LF_known",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, mean_fun = stub_known, density_fun = stub_known, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.mean_fun <- mean_fun
      private$.density_fun <- density_fun
    },
    get_mean = function(tmle_task, fold_number) {
      learner_task <- tmle_task$get_regression_task(self$name, scale = FALSE)
      preds <- self$mean_fun(learner_task)

      return(preds)
    },
    get_density = function(tmle_task, fold_number) {
      learner_task <- tmle_task$get_regression_task(self$name, scale = FALSE)
      preds <- self$density_fun(learner_task)

      outcome_type <- learner_task$outcome_type
      observed <- outcome_type$format(learner_task$Y)
      if (outcome_type$type == "binomial") {
        likelihood <- ifelse(observed == 1, preds, 1 - preds)
      } else if (outcome_type$type == "categorical") {
        unpacked <- sl3::unpack_predictions(preds)
        index_mat <- cbind(seq_along(observed), observed)
        likelihood <- unpacked[index_mat]
      } else if (outcome_type$type == "continuous") {
        likelihood <- unlist(preds)
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      return(likelihood)
    }
  ),
  active = list(
    mean_fun = function() {
      return(private$.mean_fun)
    },

    density_fun = function() {
      return(private$.density_fun)
    }
  ),
  private = list(
    .name = NULL,
    .mean_fun = NULL,
    .density_fun = NULL
  )
)
