#' Cache Likelihood values, update those values
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base
#' @importFrom assertthat assert_that is.count is.flag
#' @importFrom delayed bundle_delayed
#' @export
Likelihood_cache <- R6Class(
  classname = "Likelihood_cache",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function() {
      private$.cache <- new.env()
    },
    tasks_at_step = function(current_step) {
      self$tasks[task_uuids]
    },
    get_update_step = function(likelihood_factor, tmle_task, fold_number) {
      key <- self$key(likelihood_factor, tmle_task, fold_number)
      step_key <- sprintf("%s_%s", key, "step")
      get0(step_key, self$cache, inherits = FALSE)
    },
    key = function(likelihood_factor, tmle_task, fold_number) {
      key <- sprintf("%s_%s_%s", likelihood_factor$uuid, tmle_task$uuid, fold_number)
      return(key)
    },
    set_values = function(likelihood_factor, tmle_task, update_step = 0, fold_number, values) {
      self$cache_task(tmle_task)

      # respect likelihood factors that don't want to cache
      if (!likelihood_factor$cache) {
        return(0)
      }
      key <- self$key(likelihood_factor, tmle_task, fold_number)
      assign(key, values, self$cache)

      step_key <- sprintf("%s_%s", key, "step")
      assign(step_key, update_step, self$cache)

      return(1)
    },
    get_values = function(likelihood_factor, tmle_task, fold_number) {
      # matching_index <- self$find_match(likelihood_factor, tmle_task, fold_number)
      key <- self$key(likelihood_factor, tmle_task, fold_number)
      values <- get0(key, self$cache, inherits = FALSE)

      return(values)
    },
    cache_lf = function(likelihood_factor) {
      private$.lfs[likelihood_factor$uuid] <- likelihood_factor
    },
    cache_task = function(task) {
      if (!(task$uuid %in% names(private$.tasks))) {
        private$.tasks[[task$uuid]] <- task
      }
    }
  ),
  active = list(
    cache = function() {
      return(private$.cache)
    },
    tasks = function() {
      return(private$.tasks)
    }
  ),
  private = list(
    .tasks = list(),
    .lfs = list(),
    .cache = NULL
  )
)
