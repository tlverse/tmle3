#' Defines a TML Estimator (except for the data)
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec <- R6Class(
  classname = "tmle3_Spec",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(likelihood_override = NULL, node_types = NULL, ...) {
      private$.options <- list(likelihood_override = likelihood_override, node_types = node_types, ...)
    },
    make_tmle_task = function(data, node_list, ...) {
      setDT(data)
      
      node_types <- self$options$node_types
      if(is.null(node_types$Y)){
        Y_node <- node_list$Y
        Y_vals <- unlist(data[, Y_node, with = FALSE])
        Y_variable_type <- variable_type(x = Y_vals)  
      } else {
        Y_variable_type <- node_types$Y
      }
      
      if(!Y_variable_type$type%in%c("continuous","binomial")){
        stop("Y variable detected to be ", Y_variable_type$type,". Only continuous and binomial are supported.")
      }
      
      # bound Y if continuous
      if ((Y_variable_type$type == "continuous") && is.null(Y_variable_type$bounds)) {
        min_Y <- min(Y_vals)
        max_Y <- max(Y_vals)
        range <- max_Y - min_Y
        lower <- min_Y  #- 0.1 * range
        upper <- max_Y  #+ 0.1 * range
        Y_variable_type <- variable_type(
          type = "continuous",
          bounds = c(lower, upper)
        )
      }

      # make tmle_task
      npsem <- list(
        define_node("W", node_list$W),
        define_node("A", node_list$A, c("W"), node_types$A),
        define_node("Y", node_list$Y, c("A", "W"), Y_variable_type)
      )

      if (!is.null(node_list$id)) {
        tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
      } else {
        tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
      }

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided
      A_variable_type <- tmle_task$npsem$A$variable_type$type
      if(A_variable_type == "continuous"){
        bound <- NULL
      } else{
        bound <- 0.025
      }
      likelihood_def <- self$options$likelihood_override
      if (!is.null(likelihood_def)) {
        likelihood <- likelihood_def$train(tmle_task)
      } else {
        factor_list <- list(
          define_lf(LF_emp, "W"),
          define_lf(LF_fit, "A", learner = learner_list[["A"]], bound=bound),
          define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
        )

        likelihood_def <- Likelihood$new(factor_list)

        # fit_likelihood
        likelihood <- likelihood_def$train(tmle_task)
      }
      return(likelihood)
    },
    make_updater = function() {
      updater <- tmle3_Update$new()
    },
    make_targeted_likelihood = function(likelihood, updater) {
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
      return(targeted_likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      stop("this is a base class, try tsm_Spec_TSM_all")
    }
  ),
  active = list(
    options = function() {
      return(private$.options)
    }
  ),
  private = list(
    .options = NULL
  )
)
