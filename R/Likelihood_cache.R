# Cache Likelihood values, update those values
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
  	initialize = function(){
  		private$.cache <- data.table(
  			lf_uuid = character(),
  			task_uuid = character(),
  			update_step = integer(),
  			values = list()
  			)  			
  	},
  	column_name = function(lf_uuid, task_uuid){
  		return(sprintf("%s_%s",lf_uuid, task_uuid))
  	},
  	get_update_step = function(likelihood_factor, tmle_task){
    	req_lf_uuid <- likelihood_factor$uuid
  		req_task_uuid <- tmle_task$uuid
  		return(self$.cache[lf_uuid==req_lf_uuid & 
  									  task_uuid==req_task_uuid, 
  									  update_step])
  	},
  	get_values = function(likelihood_factor, tmle_task, update_step){
    	lf_uuid <- likelihood_factor$lf_uuid
  		task_uuid <- tmle_task$uuid
  		column_name <- self$column_name(lf_uuid, task_uuid)
  		if(column_name%in%names(self$cache)){
  			return(private$.value_cache[[column_name]])
  		} else {
  			return(NULL)
  		}
  	},
  	set_values = function(lf_uuid, task_uuid, update_step, values){
  		new_data <- list(lf_uuid=lf_uuid, 
  						 task_uuid=task_uuid, 
  						 update_step=update_step, 
  						 values=list(values))
  		
  		matching_index <- self$.cache[lf_uuid==new_data$lf_uuid & 
  									  task_uuid==new_data$task_uuid, 
  									  .I]
  		
  		if(length(matching_index)==0){
  			private$.cache <- rbindlist(list(private$.cache, new_data))
  		} else {
  			set(private$.cache, matching_index, names(private$.cache), new_data)
  		}
  		
  		return(length(matching_index))
  },
  cache_lf = function(likelihood_factor){
  	private$.lfs[likelihood_factor$uuid] <- likelihood_factor
  },
  cache_task = function(task){
  	private$.tasks[task$uuid] <- task
  }),
  active = list(
  	cache = function(){
  		return(private$.cache)
  	}
  ),
  private = list(
  	.tasks = list(),
  	.lfs = list(),
  	.cache = NULL
  )
)
  	