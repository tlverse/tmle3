#' Static Likelihood Factor
#'
#' Likelihood factor for a variable that only has one value with probability 1. This is used for static interventions.
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
#'   \code{define_lf(LF_static, name, type, value, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     \item{\code{value}}{the static value
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{value}}{the static value.}
#'     }

#'
#' @export
LF_static <- R6Class(
  classname = "LF_static",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, type = "density", value, ...) {
      super$initialize(name, ..., type = type)
      private$.value <- value
      private$.variable_type <- variable_type("constant", value)
    },
    get_mean = function(tmle_task, fold_number, ...) {
      observed <- tmle_task$get_tmle_node(self$name, include_time = T, include_id = T, expand = T, compute_risk_set = F)
      node_vars <- tmle_task$npsem[[self$name]]$variables
      #observed <- observed[, c("id", "t", node_vars), with = F]
      observed <- observed[,node_vars, with = F]

      #set(observed, , node_vars, self$value)
      pred <- unlist(observed)
      #setnames(pred, node_vars, self$name )
      return(pred)
    },
    get_density = function(tmle_task, fold_number, expand = T, ...) {
      observed <- tmle_task$get_tmle_node(self$name, include_time = T, include_id = T, expand = expand, compute_risk_set = F)
      node_vars <- tmle_task$npsem[[self$name]]$variables
      #observed <- observed[, c("id", "t", node_vars), with = F]
      observed <- observed[, node_vars, with = F]
      set(observed, , node_vars, data.table(as.numeric(self$value == unlist(observed[,node_vars, with = F ]))))
      likelihood <- unlist(observed)
      #setnames(likelihood,node_vars, self$name)
      return(likelihood)
    },
    cf_values = function(tmle_task) {
      #cf_values <- rep(self$value, tmle_task$nrow)
      node <- tmle_task$npsem[[self$name]]
      times <- node$times
      if(is.null(times)){
        times <- 0
      }
      num_times <- (length(times))
      cf_values <- rep(self$value, length(unique(tmle_task$id)) * num_times)
      return(cf_values)
    },
    sample = function(tmle_task, n_samples = NULL, fold_number = "full") {
      # TODO: fold
      # TODO: option to return task
      if (is.null(n_samples)) {
        return(tmle_task)
      }

      values <- replicate(n_samples, rep(self$value, tmle_task$nrow))
      return(values)
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
