#' TMLE fit object
#'
#' A tmle_fit object, containing initial and updated estimates, as well as data
#' about the fitting procedure. TMLE updates are calculated when the object is
#' constructed.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#'
#' @family Parameters
#'
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template tmle3_Fit_extra
#'
#' @export
#
tmle3_Fit <- R6Class(
  classname = "tmle3_Fit",
  public = list(
    initialize = function(tmle_task, likelihood, tmle_params, updater, ...) {
      if (inherits(tmle_params, "Param_base")) {
        tmle_params <- list(tmle_params)
      }
      private$.tmle_task <- tmle_task
      private$.likelihood <- likelihood
      private$.tmle_params <- tmle_params
      updater$tmle_params <- tmle_params
      private$.updater <- updater
      initial_psi <- sapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(self$tmle_task)$psi
        }
      )
      private$.initial_psi <- unlist(initial_psi)
      private$.tmle_fit(max_it)
    },
    print = function() {
      cat(sprintf("A tmle3_Fit that took %s step(s)\n", self$steps))
      print(self$summary)
    },
    set_timings = function(start_time, task_time, likelihood_time, params_time,
                           fit_time) {
      timings <- list(
        make_tmle_task = task_time - start_time,
        fit_likelihood = likelihood_time - task_time,
        define_params = params_time - likelihood_time,
        tmle_update = fit_time - params_time
      )
      private$.timings <- do.call(rbind, timings)
    }
  ),
  active = list(
    tmle_task = function() {
      return(private$.tmle_task)
    },
    likelihood = function() {
      return(private$.likelihood)
    },
    tmle_params = function() {
      return(private$.tmle_params)
    },
    tmle_param_names = function() {
      if (is.null(private$.tmle_param_names)) {
        private$.tmle_param_names <- unlist(sapply(self$tmle_params, `[[`, "name"))
      }
      return(private$.tmle_param_names)
    },
    tmle_param_types = function() {
      if (is.null(private$.tmle_param_types)) {
        private$.tmle_param_types <- sapply(self$tmle_params, `[[`, "type")
      }
      return(private$.tmle_param_types)
    },
    updater = function() {
      return(private$.updater)
    },
    steps = function() {
      return(self$updater$step_number)
    },
    ED = function() {
      ED <- private$.ED
      # names(ED) <- self$tmle_param_names
      return(ED)
    },
    initial_psi = function() {
      initial_psi <- private$.initial_psi
      # names(initial_psi) <- self$tmle_param_names
      return(initial_psi)
    },
    estimates = function() {
      estimates <- private$.estimates
      # names(estimates) <- self$tmle_param_names
      return(estimates)
    },
    summary = function() {
      return(summary_from_estimates(
        task = self$tmle_task, estimates = self$estimates,
        param_names = self$tmle_param_names,
        param_types = self$tmle_param_types,
        init_psi = self$initial_psi
      ))
    },
    timings = function() {
      return(private$.timings)
    }
  ),
  private = list(
    .tmle_task = NULL,
    .likelihood = NULL,
    .tmle_params = NULL,
    .tmle_param_names = NULL,
    .tmle_param_types = NULL,
    .updater = NULL,
    .steps = NULL,
    .ED = NULL,
    .initial_psi = NULL,
    .estimates = NULL,
    .timings = NULL,
    .tmle_fit = function(max_it = 100) {
      self$updater$update(self$likelihood, self$tmle_task)
      private$.steps <- self$updater$steps


      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(self$tmle_task, self$updater$update_fold)
        }
      )

      private$.estimates <- estimates
      private$.ED <- ED_from_estimates(estimates)
    }
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor.
#' @rdname tmle3_Fit
#' @export
#
fit_tmle3 <- tmle3_Fit$new
