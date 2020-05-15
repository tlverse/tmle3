#' Hazard Likelihood Factor Estimated using sl3.
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
    initialize = function(name, learner, is_time_variant = FALSE, ..., type = "density") {
      super$initialize(name, learner, is_time_variant = is_time_variant, ..., type = type)
    },
    get_training_task = function(tmle_task) {
       # TODO: select N_pev=A_c_prev=0 for training
      t_data <- tmle_task$get_tmle_node("t", format = TRUE)
      T_tilde_data <- tmle_task$get_tmle_node("T_tilde", format = TRUE)
      Delta_data <- tmle_task$get_tmle_node("Delta", format = TRUE)

      N_data_prev <- ifelse(t_data - 1 >= T_tilde_data & Delta_data == 1, 1, 0)
      A_c_data_prev <- ifelse(t_data - 1 >= T_tilde_data & Delta_data == 0, 1, 0)

      training_indices <- which(N_data_prev == 0 & A_c_data_prev == 0)
      # TODO: check
      tmle_task_training <- tmle_task[training_indices]
      return(tmle_task_training)
    },
    delayed_train = function(tmle_task) {
      tmle_task_training <- self$get_training_task(tmle_task)
      return(super$delayed_train(tmle_task_training))
    },
    train = function(tmle_task, learner_fit) {
      tmle_task_training <- self$get_training_task(tmle_task)
      super$train(tmle_task_training, learner_fit)
    }
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
    # get_density = function(tmle_task, fold_number) {
    #   # # transform original data into long version
    #   # short_data <-tmle_task$data
    #   # short_npsem <- tmle_task$npsem
    #   # long_data <- make_long_data(short_data, short_npsem)
    #   # long_node_list <- make_long_node_list(short_npsem)

    #   # # generate likelihood estimates for dN and dA_c
    #   # full_tmle_task <- make_long_tmle_task(long_data, long_node_list)

    #   # return(super$get_density(full_tmle_task, fold_number))

    #   learner_task <- tmle_task$get_regression_task(self$name)
    #   learner <- self$learner
    #   preds <- learner$predict_fold(learner_task, fold_number)

    #   # TODO: check
    #   predmat <- matrix(preds, nrow = tmle_task$nrow, byrow = TRUE)
    #   likelihood <- predmat
    #   return(likelihood)
    # }
  ),
  active = list(),
  private = list()
)
