#' Helper Functions for Survival Analysis
#'
#' Handles the W (covariates), A (treatment/intervention), T_tilde (time-to-event),
#' Delta (censoring indicator), t_max (the maximum time to estimate) survival data structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{survival_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname survival_tx
survival_tx_npsem <- function(node_list, variable_types = NULL) {
	# make the tmle task
	npsem <- list(
		# TODO: causal relation, handle t_max
		define_node("W", node_list$W, variable_type = variable_types$W),
		define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
    define_node("T_tilde", node_list$T_tilde, c("A", "W"), variable_type = variable_types$T_tilde),
    define_node("Delta", node_list$Delta, variable_type = variable_types$Delta),
    # TODO: check
    define_node("t", node_list$t, variable_type = variable_types$t), 
    define_node("N", node_list$N, c("A", "W", "t"), variable_type = variable_types$N),
    define_node("A_c", node_list$A_c, c("A", "W", "t"), variable_type = variable_types$A_c)   
    # define_node("dN", node_list$dN, c("A", "W", "t"), variable_type = variable_types$dN),
    # define_node("dA_c", node_list$dA_c, c("A", "W", "t"), variable_type = variable_types$dA_c)
		)

	return(npsem)
}

#' @export
#' @rdname survival_tx
survival_tx_task <- function(data, node_list, make_npsem, variable_types = NULL, ...) {
	setDT(data)

	npsem <- make_npsem(node_list, variable_types)

	if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  	} else {
    	tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  	}

  	return(tmle_task)
}

# #' @export
# #' @rdname survival_tx
# survival_tx_likelihood  <- function(tmle_task, learner_list) {
# 	# covariates
#   W_factor <- define_lf(LF_emp, "W")

#   # TODO: check if necessary
#   # treatment (bound likelihood away from 0 (and 1 if binary))
#   A_type <- tmle_task$npsem[["A"]]$variable_type
#   if (A_type$type == "continous") {
#     A_bound <- c(1 / tmle_task$nrow, Inf)
#   } else if (A_type$type %in% c("binomial", "categorical")) {
#     A_bound <- 0.025
#   } else {
#     A_bound <- NULL
#   }

#   A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

#   # outcome
#   T_tilde_factor <- define_lf(LF_fit_hazard, "T_tilde", learner = learner_list[["T_tilde"]])

#   # construct and train likelihood
#   factor_list <- list(W_factor, A_factor, T_tilde_factor)

#   likelihood_def <- Likelihood$new(factor_list)
#   likelihood <- likelihood_def$train(tmle_task)
#   return(likelihood)
# }

# #' @export
# #' @rdname survival_tx
# survival_tx_base_likelihood  <- function(tmle_task, learner_list) {
#   # covariates
#   W_factor <- define_lf(LF_emp, "W")

#   # TODO: check if necessary
#   # treatment (bound likelihood away from 0 (and 1 if binary))
#   A_type <- tmle_task$npsem[["A"]]$variable_type
#   if (A_type$type == "continous") {
#     A_bound <- c(1 / tmle_task$nrow, Inf)
#   } else if (A_type$type %in% c("binomial", "categorical")) {
#     A_bound <- 0.025
#   } else {
#     A_bound <- NULL
#   }

#   A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

#   # construct and train likelihood
#   factor_list <- list(W_factor, A_factor)

#   likelihood_def <- Likelihood$new(factor_list)
#   likelihood <- likelihood_def$train(tmle_task)
#   return(likelihood)
# }

# survival_tx_long_npsem <- function(node_list, variable_types = NULL) {
#   npsem <- list(
#     define_node("W", node_list$W, variable_type = variable_types$W),
#     define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
#     define_node("N", node_list$N, variable_type = variable_types$N),
#     define_node("A_c", node_list$A_c, variable_type = variable_types$A_c),
#     define_node("t", node_list$t, variable_type = variable_types$t),
#     # TODO: whether keep N, A_c
#     define_node("dN", node_list$dN, c("W", "A", "t"), variable_type = variable_types$dN),
#     define_node("dA_c", node_list$dA_c, c("W", "A", "t"), variable_type = variable_types$dA_c)
#     )

#   return(npsem)
# }

# make_long_tmle_task <- function(long_data, long_node_list) {
#   # TODO: function design
#   long_tmle_task <- survival_tx_task(long_data, long_node_list, survival_tx_long_npsem)

#   return(long_tmle_task)
# }

# # TODO: speed
# make_long_data <- function(short_data, short_npsem) {
#   n <- dim(short_data)[1]
#   rs <- NULL
#   t_tilde_name <- short_npsem$T_tilde$variables
#   t_max <- max(short_data[, ..t_tilde_name])

#   # change n
#   for (i in 1:n) {
#     # TODO: check
#     t_tilde <- short_data[i, ..t_tilde_name][[1]]
#     delta_name <- short_npsem$Delta$variables
#     delta <- short_data[i, ..delta_name]
#     w_name <- short_npsem$W$variables
#     w <- short_data[i, ..w_name]
#     a_name <- short_npsem$A$variables
#     a <- short_data[i, ..a_name]
#     current_df <- data.frame("t" = 1:t_max)
#     # for (t in 1:t_tilde) {
#     #   # TODO: check
#     #   current_df[t, w_name] <- w
#     #   current_df[t, a_name] <- a
#     #   # store N(t-1)
#     #   current_df[t, "N"] <- ifelse(t_tilde <= t - 1 & delta == 1, 1, 0)
#     #   current_df[t, "A_c"] <- ifelse(t_tilde <= t - 1 & delta == 0, 1, 0)

#     #   N_prev <- current_df[t, "N"]
#     #   N_t <- ifelse(t_tilde <= t & delta == 1, 1, 0)
#     #   current_df[t, "dN"] <- ifelse(N_t == 1 & N_prev == 0, 1, 0)

#     #   A_c_prev <- current_df[t, "A_c"]
#     #   A_c_t <- ifelse(t_tilde <= t & delta == 0, 1, 0)
#     #   current_df[t, "dA_c"] <- ifelse(A_c_t == 1 & A_c_prev == 0, 1, 0)
#     # }
#     current_df[, w_name] <- rep(w, t_max)
#     current_df[, a_name] <- rep(a, t_max)

#     current_df[, "N"] <- ifelse(t_tilde <= 1:t_max - 1 & rep(delta == 1, t_max), 1, 0)
#     current_df[, "A_c"] <- ifelse(t_tilde <= 1:t_max - 1 & rep(delta == 0, t_max), 1, 0)

#     N_prev <- current_df[, "N"]
#     N_t <- ifelse(t_tilde <= 1:t_max & rep(delta == 1, t_max), 1, 0)
#     current_df[, "dN"] <- ifelse(N_t == 1 & N_prev == 0, 1, 0)

#     A_c_prev <- current_df[, "A_c"]
#     A_c_t <- ifelse(t_tilde <= 1:t_max & rep(delta == 0, t_max), 1, 0)
#     current_df[, "dA_c"] <- ifelse(A_c_t == 1 & A_c_prev == 0, 1, 0)

#     if (i == 1) {
#       rs <- current_df
#     } else {
#       rs <- rbind(rs, current_df)
#     }
#   }
#   return(rs)
# }

# make_long_node_list <- function(short_npsem) {
#   w_name <- short_npsem$W$variables
#   a_name <- short_npsem$A$variables
#   node_list <- list(W = w_name, A = a_name, N = "N", A_c = "A_c", t = "t", dN = "dN", dA_c = "dA_c")
#   return(node_list)
# }

#' @export
#' @rdname survival_tx
survival_tx_likelihood  <- function(tmle_task, learner_list) {
  # covariates
  # TODO: whether remove duplicates for LF_emp
  W_factor <- define_lf(LF_emp, "W")

  # TODO: check if necessary
  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

  # npsem <- tmle_task$npsem
  # delta_name <- npsem$Delta$variables
  # # TODO: check
  # Delta <- tmle_task$data[[delta_name]]
  # t_tilde_name <- npsem$T_tilde$variables
  # T_tilde <- tmle_task$data[[t_tilde_name]]

  # dN_learner <- learner_list[["dN"]]
  # # TODO: check
  # dN_factor <- define_lf(LF_fit_hazards, "dN", 
  #   learner = make_learner(Lrnr_conditional_hazards, "dN", Delta, T_tilde, dN_learner))
  # dA_c_learner <- learner_list[["dA_c"]]
  # dA_c_factor <- define_lf(LF_fit_hazards, "dA_c", 
  #   learner = make_learner(Lrnr_conditional_hazards, "dA_c", Delta, T_tilde, dA_c_learner))

  # TODO: modify get_regression_task and LF_fit for time variance
  N_factor <- define_lf(LF_fit_hazards, "N", learner = learner_list[["N"]], is_time_variant = TRUE)
  A_c_factor <- define_lf(LF_fit_hazards, "A_c", learner = learner_list[["A_c"]], is_time_variant = TRUE)

  factor_list <- list(W_factor, A_factor, N_factor, A_c_factor)

  # TODO: check
  likelihood_def <- Likelihood$new(factor_list)
  # # TODO: select N=A_c=0 for training, would this affect gA and gW? only gW
  # N_data <- tmle_task$get_tmle_node("N", format = TRUE)
  # A_c_data <- tmle_task$get_tmle_node("A_c", format = TRUE)
  # training_indices <- which(N_data == 0 & A_c_data == 0)
  # # training_data <- tmle_task$data[training_indices,]
  # # TODO: check
  # # tmle_task_training <- tmle_task$next_in_chain(data = training_data)
  # tmle_task_training <- tmle_task[training_indices]
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

# # TODO: optimize
# convert_node_name <- function(node, t_max) {
#   rs <- unlist(lapply(seq(t_max), function(i) {
#       paste(node, as.character(i), sep = "_")
#       }))
#   return(rs)
# }
