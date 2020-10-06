







sample_all <- function(observed_task, tlik, num_samples = 5000, baseline_var = "W", time_ordering) {
  mc_data <- as.data.table(tlik$factor_list[[baseline_var]]$sample(, num_samples))
  set(mc_data, , setdiff(names(observed_task$data), "W"), (0))
  mc_data$t <- 0
  mc_data$id <- 1:nrow(mc_data)
  mc_task <- tmle3_Task$new(mc_data, observed_task$npsem, id = "id", t = "t", long_format  = observed_task$long_format)
  time_ordering <- setdiff(time_ordering, baseline_var)
  for(node in time_ordering) {
    mc_task <- sample_from_node(mc_task, tlik, node, observed_task)
  }
  return(mc_task)
}


# Function that samples from a single node of a targeted likelihood
# assumes that mc_task has all the columns of the original data. Rows will be added if this nodes time is not yet in data.
sample_from_node <- function(mc_task, tlik, node, observed_task) {
  data <- mc_task$data
  var <- observed_task$npsem[[node]]$variables
  time <- observed_task$npsem[[node]]$time
  levels <- sort(unique(observed_task$get_tmle_node(node)))

  if(!(time %in% unique(data$t))){
    data1 <- data[data$t==0]
    data2<- copy(data1)
    set(data2, , setdiff(colnames(data2), "id"),  0)
    data2$t <- time
    data <- rbind(data, data2)

  }
  expand_data <- function(data, levels, time, var){
    expanded_data <- rbindlist(lapply(levels, function(level) {
      data <- copy(data)
      data$level <- level
      data$trueid <- data$id
      set(data, which(data$t == time), var, level)

      return(data)
    }))
    expanded_data$id <- paste0(expanded_data$trueid, "_", expanded_data$level)
    setkey(expanded_data, id, t)
    return(expanded_data)
  }

  expanded_data <- expand_data(data, levels, time, var)
  expanded_task <- tmle3_Task$new(expanded_data, observed_task$npsem, id = "id", t = "t", long_format  = observed_task$long_format, summary_measure_columns = c("trueid", observed_task$summary_measure_columns))

  setattr(expanded_task, "target_nodes", c(node))
  if(inherits(tlik, "Targeted_Likelihood")) {
    tlik$sync_task(expanded_task, check = F)
  }


  node_liks <- data.table(trueid = expanded_task$data$trueid, expanded_task$get_tmle_node(node, include_id = T)
                          , lik  = as.vector(tlik$get_likelihood(expanded_task, node) ))

  #wide_liks <- (dcast(node_liks, as.formula(paste0("trueid ~ ", var)), , value.var = "lik"))


  #sampled <- (as.data.table(t(apply(wide_liks, 1, function(v, levels){
   # c(v[1], sample(levels, 1, prob = v[-1]))
  #}, levels = levels))))
  #setnames(sampled, c("id", node))

  sampled <- node_liks[, sample(.SD[[1]], 1, prob = lik) , by = "trueid", .SDcols = var]
  setnames(sampled, c("id", node))
  setkey(sampled, id)
  sampled$t <- time
  sampled$id <- as.factor(sampled$id)

  mc_task_new <- mc_task$generate_counterfactual_task(UUIDgenerate(), sampled)
  return(mc_task_new)
}






























#A general sampler from likelihood objects
#' @importFrom AR AR.Sim
#' @param likelihood The likelihood object to sample from
#' @param time_ordering A character vector of node names which represents the time ordering to sample form.
#' @param use_lf A character vector of node names whose likelihood factor should be used for sampling (e.g. nodes estimated with LF_emp)
#' @export

Sampler <- R6Class(
  classname = "Sampler",
  portable = TRUE,
  class = TRUE,
  active = list(
  params = function(){
    private$.params
  },
  time_ordering = function(){
    self$params$time_ordering
  },
  likelihood = function(){
    self$params$likelihood
  }
  ),
  public = list(
    initialize = function(likelihood, time_ordering, use_lf = c()){
      params <- sl3::args_to_list()
      private$.params <- params
    },
    sample = function(tmle_task, start_node, end_node){
      start_index <- which(self$time_ordering == start_node)
      end_index <- which(self$time_ordering == end_node)
      if(end_index < start_index){
        stop("Start and end node do not satisfy time ordering.")
      }
      for(i in start_index:end_index){
        node <- self$time_ordering[[i]]
        tmle_task <- self$next_in_chain(tmle_task, node)
      }
      return(tmle_task)

    },
    compute_conditional_mean = function(tmle_task, start_node, end_node, num_iter = 100){
      # Computes conditional mean of end_node given the (strict) past of start_node via monte carlo simulation
      # The sampling begins with start_node.
      outcomes = matrix(nrow = length(unique(tmle_task$id)), ncol = num_iter)
      for(i in 1:num_iter){
        new_task <- self$sample(tmle_task, start_node, end_node)
        outcomes[,i] <- new_task$get_tmle_node(end_node)[,end_node,with=F][[1]]
      }
      return(rowMeans(outcomes))
    },
    next_in_chain = function(tmle_task, node){
      #Takes a task and node name and then generated a new task
      # with the node values replaced with sampled values.
      samples <- data.table(self$sample_from_node(tmle_task, node, 1))
      setnames(samples, node)
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), samples)
      return(cf_task)
    },
    sample_from_node = function(tmle_task, node, n_samples = 1, use_LF_factor = node %in% self$params$use_lf, fold_number = "full"){
      #Samples from a node. Returns matrix of n by n_samples of values.

      if(use_LF_factor){
        return(self$likelihood$factor_list[[node]]$sample(tmle_task, n_samples, fold_number))
      }
      outcome_type <- tmle_task$npsem[[node]]$variable_type
      times_to_pool <- tmle_task$npsem[[node]]$times_to_pool

      num_id <- length(unique(tmle_task$id))

      num_rows <- ifelse(is.null(times_to_pool), num_id, length(times_to_pool)*num_id)
      if(outcome_type$type == "binomial"){
        cf_outcome <- data.table(A=rep(1, num_rows))

        setnames(cf_outcome, node)
        cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = cf_outcome)
        p <- self$likelihood$get_likelihood(cf_task, node)[, node, with = F][[1]]

        values <- matrix(unlist(lapply(1:num_rows, function(i) {
          as.vector(rbinom(n_samples, 1, p[i]))
        })), ncol = num_rows)
      }
      else if (outcome_type$type == "categorical"){
        levels <- outcome_type$levels
        cf_tasks <- lapply(levels, function(level){
          cf_outcome <- data.table(rep(level,num_rows))
          setnames(cf_outcome, node)
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_outcome)
        })
        probs <- as.matrix(lapply(cf_tasks, function(cf_task){
          p <- self$likelihood$get_likelihood(cf_task, node)[, node, with = F][[1]]
        }))
        values <- apply(probs, 1, function(p){
            apply(
              rmultinom(n_samples, 1, p) == 1, 2,
              function(onehots) levels[which(onehots)]
            )
          })
      }
      else if (outcome_type$type == "continuous") {
        if(!is.null(times_to_pool)){
          stop("Sampling from continuous variables is not supported for pooled time.")
        }
        outcome <- tmle_task$get_tmle_node(node)
        outcome_stripped <- outcome[, node, with = F][[1]]
        ids <- outcome$id
        values <- matrix(nrow = n_samples, ncol = num_rows)
        for (id in ids) {
          subject <- tmle_task[which(tmle_task$id == id)]
          f_X <- function(a) {

            #TODO might be more efficient to pass the regression task to likelihood so we dont recompute
            cf_data <- data.table(a)
            print(cf_data)
            setnames(cf_data, names(cf_data), node)
            subject_a <- subject$generate_counterfactual_task(UUIDgenerate(), cf_data)
            likelihood <- self$likelihood$get_likelihood(subject_a, node)[, node, with = F][[1]]
            return(likelihood)
          }
          samples <- AR::AR.Sim(n_samples, f_X, Y.dist = "norm", Y.dist.par = c(mean(outcome_stripped), var(outcome_stripped)),
                            xlim = c(min(outcome_stripped), max(outcome_stripped))
          )
          values[, i] <- samples
        }

      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      values <- t(values)
      return(values)

    }
  ),
  private = list(
    .params = NULL
  )
)
