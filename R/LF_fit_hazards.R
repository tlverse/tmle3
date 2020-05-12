#' Hazard Likelihood Factor Estimated from Transformed Data using sl3.
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
LF_fit_hazards <- R6Class(
  classname = "LF_fit_hazards",
  portable = TRUE,
  class = TRUE,
  inherit = LF_fit,
  public = list(
    initialize = function(name, learner, ..., type = "density") {
      super$initialize(name, learner, ..., type = type)
    },
    # delayed_train = function(tmle_task) {
    #   # just return prefit learner if that's what we have
    #   # otherwise, make a delayed fit and return that
    #   if (self$learner$is_trained) {
    #     return(self$learner)
    #   }

    #   # transform original data into long version
    #   short_data <-tmle_task$data
    #   short_npsem <- tmle_task$npsem
    #   long_data <- make_long_data(short_data, short_npsem)
    #   long_node_list <- make_long_node_list(short_npsem)

    #   # make long tmle task
    #   selected_long_data <- long_data[which(long_data["N"] == 0 & long_data["A_c"] == 0),]
    #   # TODO: variable type
    #   long_tmle_task <- make_long_tmle_task(selected_long_data, long_node_list)

    #   return(super$delayed_train(long_tmle_task))
    # },
    # train = function(tmle_task, learner_fit) {
    #   # transform original data into long version
    #   short_data <-tmle_task$data
    #   short_npsem <- tmle_task$npsem
    #   long_data <- make_long_data(short_data, short_npsem)
    #   long_node_list <- make_long_node_list(short_npsem)

    #   # make long tmle task
    #   selected_long_data <- long_data[which(long_data["N"] == 0 & long_data["A_c"] == 0),]
    #   # TODO: variable type
    #   long_tmle_task <- make_long_tmle_task(selected_long_data, long_node_list)

    #   super$train(long_tmle_task, learner_fit)
    # },
    get_density = function(tmle_task, fold_number) {
      # # transform original data into long version
      # short_data <-tmle_task$data
      # short_npsem <- tmle_task$npsem
      # long_data <- make_long_data(short_data, short_npsem)
      # long_node_list <- make_long_node_list(short_npsem)

      # # generate likelihood estimates for dN and dA_c
      # full_tmle_task <- make_long_tmle_task(long_data, long_node_list)

      # return(super$get_density(full_tmle_task, fold_number))

      learner_task <- tmle_task$get_regression_task(self$name)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)

      # TODO: check
      predmat <- matrix(preds, nrow = tmle_task$nrow, byrow = TRUE)
      likelihood <- predmat
      return(likelihood)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL
  )
)
