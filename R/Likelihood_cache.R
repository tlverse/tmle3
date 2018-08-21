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
      private$.cache <- data.table(
        lf_uuid = character(),
        task_uuid = character(),
        update_step = integer(),
        cv_fold = integer(),
        values = list()
      )
    },
    find_match = function(likelihood_factor, tmle_task, search_cv_fold) {
      self$cache_task(tmle_task)
      matching_index <- private$.cache[, which(lf_uuid == likelihood_factor$uuid &
        task_uuid == tmle_task$uuid & cv_fold == search_cv_fold)]
      return(matching_index)
    },
    tasks_at_step = function(current_step) {
      matching_index <- private$.cache[, which(update_step == current_step)]
      task_uuids <- unique(private$.cache$task_uuid[matching_index])
      self$tasks[task_uuids]
    },
    get_update_step = function(likelihood_factor, tmle_task, cv_fold) {
      matching_index <- self$find_match(likelihood_factor, tmle_task, cv_fold)
      if (length(matching_index) == 0) {
        return(NULL)
      } else {
        return(private$.cache$update_step[[matching_index]])
      }
    },
    get_values = function(likelihood_factor, tmle_task, cv_fold) {
      matching_index <- self$find_match(likelihood_factor, tmle_task, cv_fold)

      if (length(matching_index) == 0) {
        return(NULL)
      } else {
        return(private$.cache$values[[matching_index]])
      }
    },
    set_values = function(likelihood_factor, tmle_task, update_step = 0, cv_fold, values) {
      new_data <- list(
        lf_uuid = likelihood_factor$uuid,
        task_uuid = tmle_task$uuid,
        update_step = update_step,
        cv_fold = cv_fold,
        values = list(values)
      )

      matching_index <- self$find_match(likelihood_factor, tmle_task, cv_fold)

      if (length(matching_index) == 0) {
        private$.cache <- rbindlist(list(private$.cache, new_data))
      } else {
        set(
          private$.cache,
          matching_index,
          names(private$.cache),
          new_data
        )
      }

      return(length(matching_index))
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
