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
        values = list()
      )
    },
    column_name = function(lf_uuid, task_uuid) {
      return(sprintf("%s_%s", lf_uuid, task_uuid))
    },
    find_match = function(search_lf_uuid, search_task_uuid){
      matching_index <- private$.cache[,which(lf_uuid == search_lf_uuid &
                                                task_uuid == search_task_uuid)] 
      return(matching_index)
    },
    get_update_step = function(lf_uuid, task_uuid) {
      matching_index <- self$find_match(lf_uuid,task_uuid)
      if (length(matching_index) == 0) {
        return(NULL)
      } else {
        return(private$.cache$update_step[[matching_index]])  
      }
    },
    get_values = function(lf_uuid, task_uuid, update_step = 0) {
      # 	lf_uuid <- likelihood_factor$uuid
      # task_uuid <- tmle_task$uuid
      # column_name <- self$column_name(lf_uuid, task_uuid)
      # if(column_name%in%names(self$cache)){
      # 	return(private$.value_cache[[column_name]])
      # } else {
      # 	return(NULL)
      # }

      matching_index <- self$find_match(lf_uuid, task_uuid)
      
      if (length(matching_index) == 0) {
        return(NULL)
      } else {
        return(private$.cache$values[[matching_index]])
      }
      
    },
    set_values = function(lf_uuid, task_uuid, update_step = 0, values) {
      new_data <- list(
        lf_uuid = lf_uuid,
        task_uuid = task_uuid,
        update_step = update_step,
        values = list(values)
      )
      
      matching_index <- self$find_match(lf_uuid, task_uuid)
      
      if (length(matching_index) == 0) {
        private$.cache <- rbindlist(list(private$.cache, new_data))
      } else {
        set(private$.cache,
            matching_index,
            names(private$.cache),
            new_data)
      }
      
      return(length(matching_index))
    },
    cache_lf = function(likelihood_factor) {
      private$.lfs[likelihood_factor$uuid] <- likelihood_factor
    },
    cache_task = function(task) {
      private$.tasks[task$uuid] <- task
    }
  ),
  active = list(
    cache = function() {
      return(private$.cache)
    }
  ),
  private = list(
    .tasks = list(),
    .lfs = list(),
    .cache = NULL
  )
)