#' Defines an update procedure (submodel+loss function) for survival data
#'
#' Current Limitations:
#' loss function and submodel are hard-coded (need to accept arguments for these)
#' @section Constructor:
#'   \code{define_param(maxit, cvtmle, one_dimensional, constrain_step, delta_epsilon, verbose)}
#'
#'   \describe{
#'     \item{\code{maxit}}{The maximum number of update iterations
#'     }
#'     \item{\code{cvtmle}}{If \code{TRUE}, use CV-likelihood values when
#'        calculating updates.
#'     }
#'     \item{\code{one_dimensional}}{If \code{TRUE}, collapse clever covariates
#'        into a one-dimensional clever covariate scaled by the mean of their
#'        EIFs.
#'     }
#'     \item{\code{constrain_step}}{If \code{TRUE}, step size is at most
#'        \code{delta_epsilon} (it can be smaller if a smaller step decreases
#'        the loss more).
#'     }
#'     \item{\code{delta_epsilon}}{The maximum step size allowed if
#'        \code{constrain_step} is \code{TRUE}.
#'     }
#'     \item{\code{convergence_type}}{The convergence criterion to use: (1)
#'        \code{"scaled_var"} corresponds to sqrt(Var(D)/n)/logn (the default)
#'        while (2) \code{"sample_size"} corresponds to 1/n.
#'     }
#'     \item{\code{fluctuation_type}}{Whether to include the auxiliary covariate
#'        for the fluctuation model as a covariate or to treat it as a weight.
#'        Note that the option \code{"weighted"} is incompatible with a
#'        multi-epsilon submodel (\code{one_dimensional = FALSE}).
#'     }
#'     \item{\code{verbose}}{If \code{TRUE}, diagnostic output is generated
#'        about the updating procedure.
#'     }
#'     }
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Update_survival <- R6Class(
  classname = "tmle3_Update_survival",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Update,
  public = list(
    initialize = function(maxit = 100, cvtmle = TRUE, one_dimensional = FALSE,
                              constrain_step = FALSE, delta_epsilon = 1e-4,
                              convergence_type = c("scaled_var", "sample_size"),
                              fluctuation_type = c("standard", "weighted"),
                              use_best = FALSE,
                              verbose = FALSE,
                              fit_method = "l2",
                              clipping = 1e0) {
      super$initialize(maxit, cvtmle, one_dimensional,
                              constrain_step , delta_epsilon,
                              convergence_type,
                              fluctuation_type,
                              use_best,
                              verbose)
     private$.fit_method <- fit_method
     private$.clipping <- clipping
    },
    # collapse_covariates = function(estimates, clever_covariates) {
    #   ED <- ED_from_estimates(estimates)
    #   EDnormed <- ED / norm(ED, type = "2")
    #   collapsed_covariate <- clever_covariates %*% EDnormed

    #   return(collapsed_covariate)
    # },
    # update_step = function(likelihood, tmle_task, fold_number = "full") {

    #   # get new submodel fit
    #   all_submodels <- self$generate_submodel_data(
    #     likelihood, tmle_task,
    #     fold_number
    #   )
    #   new_epsilons <- self$fit_submodels(all_submodels)

    #   # update likelihoods
    #   likelihood$update(new_epsilons, self$step_number, fold_number)

    #   if (fold_number != "full") {
    #     # update full fit likelihoods if we haven't already
    #     likelihood$update(new_epsilons, self$step_number, "full")
    #   }
    #   # increment step count
    #   private$.step_number <- private$.step_number + 1
    # },
    # generate_submodel_data = function(likelihood, tmle_task,
    #                                       fold_number = "full") {
    #   update_nodes <- self$update_nodes

    #   # TODO: support not getting observed for case where we're applying
    #   #       updates instead of fitting them
    #   clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
    #     tmle_param$clever_covariates(tmle_task, fold_number)
    #   })

    #   observed_values <- lapply(update_nodes, tmle_task$get_tmle_node)

    #   all_submodels <- lapply(update_nodes, function(update_node) {
    #     node_covariates <- lapply(clever_covariates, `[[`, update_node)
    #     covariates_dt <- do.call(cbind, node_covariates)
    #     if (self$one_dimensional) {
    #       observed_task <- likelihood$training_task
    #       estimates <- lapply(self$tmle_params, function(tmle_param) {
    #         tmle_param$estimates(observed_task, fold_number)
    #       })
    #       covariates_dt <- self$collapse_covariates(estimates, covariates_dt)
    #     }

    #     observed <- tmle_task$get_tmle_node(update_node)
    #     initial <- likelihood$get_likelihood(
    #       tmle_task, update_node,
    #       fold_number
    #     )

    #     # scale observed and predicted values for bounded continuous
    #     observed <- tmle_task$scale(observed, update_node)
    #     initial <- tmle_task$scale(initial, update_node)

    #     # protect against qlogis(1)=Inf
    #     initial <- bound(initial, 0.005)

    #     submodel_data <- list(
    #       observed = observed,
    #       H = covariates_dt,
    #       initial = initial
    #     )
    #   })

    #   names(all_submodels) <- update_nodes

    #   return(all_submodels)
    # },
    norm_l2 = function(beta) {
      return(sqrt(sum(beta^2)))
    },
    # TODO: check
    fit_submodel = function(submodel_data) {
      if (self$fit_method == "l2") {
        # TODO: check
        # print("l2")
        # mean_eic <- self$get_mean_eic(self$update_fold)
        epsilon_n <- tryCatch({
          alpha <- 0; norm_func <- self$norm_l2; lambda.min.ratio = 1e-2

          ind <- 1
          while (ind == 1) {
            submodel_fit <- glmnet::glmnet(
              x = submodel_data$H,
              y = submodel_data$observed,
              offset = qlogis(submodel_data$initial),
              family = "binomial",
              alpha = alpha,
              standardize = FALSE,
              intercept = FALSE,
              lambda.min.ratio = lambda.min.ratio,
              # nlambda = 2e2
              nlambda = 1e2
              # TODO: check
              # penalty.factor = 1/abs(mean_eic)
              )
            norms <- apply(submodel_fit$beta, 2, norm_func)
            ind <- max(which(norms <= self$clipping))
            if (ind > 1) break
            
            fit_lambda <- submodel_fit$lambda
            
            if(fit_lambda==1){
              stop("only one lambda could be fit")
            }
            
            # try to estimate what the correct lambda value is and go a bit beyond that
            norm_ratio <- self$clipping/norms[2]
            lambda_guess <- fit_lambda[1]-norm_ratio*(fit_lambda[1]-fit_lambda[2])
            lambda_min_ratio <- 0.8*lambda_guess/fit_lambda[1]
            # lambda.min.ratio <- sort(submodel_fit$lambda, decreasing = TRUE)[2] / max(submodel_fit$lambda)
          }
          epsilon_n <- submodel_fit$beta[, ind]
        }, error = function(e) {
          # TODO: check
          print(e)
          return(rep(0, ncol(submodel_data$H)))
        })
        epsilon <- epsilon_n
        # TODO: check if necessary
        # NOTE: this protects against collinear covariates
        # (which we don't care about, we just want an update)
        epsilon[is.na(epsilon)] <- 0

        if (self$verbose) {
          cat(sprintf("epsilon: %e ", epsilon))
        }
      } else {
        # TODO: check
        # print("classic")
        epsilon <- super$fit_submodel(submodel_data)
      } 
    
      return(epsilon)
    },
    # TODO: check
    get_mean_eic = function(update_fold) {
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number = update_fold)
        }
      )
      IC <- do.call(cbind, lapply(estimates, `[[`, "IC"))
      mean_eic <- colMeans(IC)
      return(mean_eic)
    },
    get_mean_eic_inner_prod = function(update_fold, tmle_task) {
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number = update_fold)
        }
      )
      IC <- do.call(cbind, lapply(estimates, `[[`, "IC"))
      mean_eic <- colMeans(IC)
      return(abs(sqrt(sum(mean_eic ^ 2))))
    }
    # update = function(likelihood, tmle_task) {
    #   update_fold <- self$update_fold
    #   maxit <- private$.maxit
    #   # TODO: check
    #   eic_list <- c(self$get_mean_eic_inner_prod(update_fold, tmle_task))
    #   count <- 0
    #   tol <- 0
    #   for (steps in seq_len(maxit)) {
    #     self$update_step(likelihood, tmle_task, update_fold)
    #     # if (self$check_convergence(tmle_task, update_fold)) {
    #     #   break
    #     # }
    #     # TODO: check
    #     mean_eic_inner_prod_current <- self$get_mean_eic_inner_prod(update_fold, tmle_task)
    #     if (mean_eic_inner_prod_current >= tail(eic_list, n=1)) {
    #       count = count + 1
    #       if (count > tol) {
    #         break
    #       }
    #     }
    #     eic_list <- c(eic_list, mean_eic_inner_prod_current)
    #   }
    #   # TODO: check
    #   return(eic_list)
    # }
    # fit_submodels = function(all_submodels) {
    #   all_epsilon <- lapply(all_submodels, self$fit_submodel)

    #   names(all_epsilon) <- names(all_submodels)
    #   private$.epsilons <- c(private$.epsilons, list(all_epsilon))

    #   return(all_epsilon)
    # },
    # submodel = function(epsilon, initial, H) {
    #   plogis(qlogis(initial) + H %*% epsilon)
    # },
    # loss_function = function(estimate, observed) {
    #   -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
    # },
    # apply_submodel = function(submodel_data, epsilon) {
    #   self$submodel(epsilon, submodel_data$initial, submodel_data$H)
    # },
    # apply_update = function(tmle_task, likelihood, fold_number, all_epsilon) {
    #   update_nodes <- self$update_nodes

    #   # get submodel data for all nodes
    #   all_submodel_data <- self$generate_submodel_data(
    #     likelihood, tmle_task,
    #     fold_number
    #   )

    #   # apply update to all nodes
    #   updated_likelihoods <- lapply(update_nodes, function(update_node) {
    #     submodel_data <- all_submodel_data[[update_node]]
    #     epsilon <- all_epsilon[[update_node]]
    #     updated_likelihood <- self$apply_submodel(submodel_data, epsilon)

    #     # un-scale to handle bounded continuous
    #     updated_likelihood <- tmle_task$unscale(
    #       updated_likelihood,
    #       update_node
    #     )
    #   })
    #   names(updated_likelihoods) <- update_nodes

    #   return(updated_likelihoods)
    # },
    # TODO: check
    # check_convergence = function(tmle_task, fold_number = "full") {
    #   estimates <- lapply(
    #     self$tmle_params,
    #     function(tmle_param) {
    #       tmle_param$estimates(tmle_task, fold_number = fold_number)
    #     }
    #   )

    #   if (self$convergence_type == "scaled_var") {
    #     # NOTE: the point of this criterion is to avoid targeting in an overly
    #     #       aggressive manner, as we simply need check that the following
    #     #       condition is met |P_n D*| / SE(D*) =< max(1/log(n), 1/10)
    #     IC <- do.call(cbind, lapply(estimates, `[[`, "IC"))
    #     se_Dstar <- sqrt(apply(IC, 2, var) / tmle_task$nrow)
    #     ED_threshold <- se_Dstar / min(log(tmle_task$nrow), 10)
    #   } else if (self$convergence_type == "sample_size") {
    #     ED_threshold <- 1 / tmle_task$nrow
    #   }

    #   # get |P_n D*| of any number of parameter estimates
    #   ED <- ED_from_estimates(estimates)
    #   ED_criterion <- abs(ED)

    #   if (self$verbose) {
    #     cat(sprintf("max(abs(ED)): %e\n", ED_criterion))
    #   }
    #   return(all(ED_criterion <= ED_threshold))
    # },
    # update = function(likelihood, tmle_task) {
    #   update_fold <- self$update_fold
    #   maxit <- private$.maxit
    #   for (steps in seq_len(maxit)) {
    #     self$update_step(likelihood, tmle_task, update_fold)
    #     if (self$check_convergence(tmle_task, update_fold)) {
    #       break
    #     }
    #   }
    # },
    # register_param = function(new_params) {
    #   if (inherits(new_params, "Param_base")) {
    #     new_params <- list(new_params)
    #   }
    #   private$.tmle_params <- c(private$.tmle_params, new_params)
    #   new_update_nodes <- unlist(lapply(new_params, `[[`, "update_nodes"))
    #   private$.update_nodes <- unique(c(
    #     private$.update_nodes,
    #     new_update_nodes
    #   ))
    # }
  ),
  active = list(
    fit_method = function() {
      return(private$.fit_method)
    },
    clipping = function() {
      return(private$.clipping)
    }
  ),
  private = list(
    .fit_method = NULL,
    .clipping = NULL
  )
)
