#' @export
resample <- function(likelihood, fun, n_samples = NULL, bootstrap=FALSE, convergence = "sample_size"){
	original_n <- likelihood$training_task$nrow
	if(is.null(n_samples)){
		n_samples <- 100
		# TODO: add support for convergence criteria
	}
	
	# TODO: parallelize using future or dopar
	if (!bootstrap) {
	  all_results <- foreach(i=1:original_n)%do%{
	    tmle_task <- likelihood$training_task[i]
	    sample_task <- likelihood$sample(1*n_samples, resample_marginal=bootstrap)
	    result <- fun(sample_task)
	    return(result)
	  }
	} else {
	  # TODO: when bootstrap=TRUE
	  all_results <- NULL
	}

	return(all_results)
}
