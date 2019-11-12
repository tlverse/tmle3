#' @export
resample <- function(likelihood, fun, n_samples = NULL, bootstrap=FALSE, convergence = "sample_size"){
	original_n <- likelihood$training_task$nrow
	if(is.null(n_samples)){
		n_samples <- 100
		# TODO: add support for convergence criteria
	}
	
	# TODO: parallelize using future or dopar
	all_results <- foreach(i=1:n_samples)%do%{
		tmle_task <- likelihood$sample(original_n, resample_marginal=bootstrap)
		result <- fun(tmle_task)
		return(result)
	}

	return(all_results)
}
