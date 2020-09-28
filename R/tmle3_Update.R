#' Defines an update procedure (submodel+loss function)
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
#'     \item{\code{use_best}}{If \code{TRUE}, the final updated likelihood is set to the
#'        likelihood that minimizes the ED instead of the likelihood at the last update
#'        step.
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
tmle3_Update <- R6Class(
  classname = "tmle3_Update",
  portable = TRUE,
  class = TRUE,
  public = list(
    # TODO: change maxit for test
    initialize = function(maxit = 100, cvtmle = TRUE, one_dimensional = FALSE,
                          constrain_step = FALSE, delta_epsilon = 1e-4,
                          convergence_type = c("scaled_var", "sample_size"),
                          fluctuation_type = c("standard", "weighted"),
                          optim_delta_epsilon = TRUE,
                          use_best = FALSE,
                          verbose = FALSE) {
      private$.maxit <- maxit
      private$.cvtmle <- cvtmle
      private$.one_dimensional <- one_dimensional
      private$.constrain_step <- constrain_step
      private$.delta_epsilon <- delta_epsilon
      private$.convergence_type <- match.arg(convergence_type)
      private$.fluctuation_type <- match.arg(fluctuation_type)
      private$.optim_delta_epsilon <- optim_delta_epsilon
      private$.use_best <- use_best
      private$.verbose <- verbose
    },
    collapse_covariates = function(estimates, clever_covariates) {

      ED <- ED_from_estimates(estimates)

      EDnormed <- ED / norm(ED, type = "2")
      collapsed_covariate <- clever_covariates %*% EDnormed

      return(collapsed_covariate)
    },
    update_step = function(likelihood, tmle_task, fold_number = "full") {
      update_nodes <- self$update_nodes

      current_step <- self$step_number + 1

      # initialize epsilons for this step
      na_epsilons <- as.list(rep(NA, length(update_nodes)))
      names(na_epsilons) <- update_nodes
      private$.epsilons[[current_step]] <- na_epsilons

      for (update_node in update_nodes) {
        # get new submodel fit
        submodel_data <- self$generate_submodel_data(
          likelihood, tmle_task,
          fold_number, update_node,
          drop_censored = TRUE,
          for_fitting = T
        )

        new_epsilon <- self$fit_submodel(submodel_data)

        # update likelihoods
        likelihood$update(new_epsilon, current_step, fold_number, update_node)

        if (fold_number != "full") {
          # update full fit likelihoods if we haven't already
          likelihood$update(new_epsilon, self$step_number, "full", update_node)
        }

        private$.epsilons[[current_step]][[update_node]] <- new_epsilon
      }

      # update step number
      private$.step_number <- current_step
    },
    generate_submodel_data = function(likelihood, tmle_task,
                                      fold_number = "full",
                                      update_node = "Y",
                                      drop_censored = FALSE, for_fitting = F) {


      if(!(inherits(likelihood, "Targeted_Likelihood"))) {
        submodel_type <- "logistic"
      } else {
        submodel_type <- likelihood$submodel_type(update_node)
      }

      submodel_info <- submodel_spec(submodel_type)
      # TODO: change clever covariates to allow only calculating some nodes
      clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
        # Assert that it supports the submodel type
        tmle_param$supports_submodel_type(submodel_type, update_node)
        #formal_args <- names(formals(tmle_param$clever_covariates))

        # For backwards compatibility:
        # In future, clever covariate functions should accept a "node" and "submodel_type" argument.
        args <- list(for_fitting = for_fitting, submodel_type = submodel_type, fold_number = fold_number, tmle_task = tmle_task
             ,node = update_node)
        return(sl3:::call_with_args(tmle_param$clever_covariates, args))
        # if("for_fitting" %in% formal_args) {
        #   return(tmle_param$clever_covariates(tmle_task, fold_number, for_fitting = for_fitting))
        # }
        # else if(all(c("submodel_type", "node") %in% formal_args)){
        #   return(tmle_param$clever_covariates(tmle_task, fold_number, submodel_type = submodel_type, node = update_node))
        # }
        # else if("submodel_type" %in% formal_args){
        #   return(tmle_param$clever_covariates(tmle_task, fold_number, submodel_type = submodel_typee))
        # }
        # else if("node" %in% formal_args){
        #   return(tmle_param$clever_covariates(tmle_task, fold_number, node = update_node))
        # }
        #  else {
        #   return(tmle_param$clever_covariates(tmle_task, fold_number))
        # }
      })

      node_covariates <- lapply(clever_covariates, `[[`, update_node)
      # Get EDs if present. Only for training task
      if(self$one_dimensional & for_fitting) {
        IC <- lapply(clever_covariates, `[[`, "IC")
        IC <- do.call(cbind, lapply(IC, `[[`, update_node) )
        if(is.null(IC)) {
          IC <- lapply(private$.current_estimates, `[[`, "IC")
          IC <- do.call(cbind, IC)
        }
      }

      covariates_dt <- do.call(cbind, node_covariates)

      # if (self$one_dimensional) {
      #   #TODO this should happen in fit_submodel so we can store epsiln
      #   observed_task <- likelihood$training_task
      #   covariates_dt <- self$collapse_covariates(self$current_estimates, covariates_dt)
      #
      #   }

      observed <- tmle_task$get_tmle_node(update_node)
      initial <- likelihood$get_likelihood(
        tmle_task, update_node,
        fold_number)

      # scale observed and predicted values for bounded continuous
      observed <- tmle_task$scale(observed, update_node)
      initial <- tmle_task$scale(initial, update_node)


      # protect against qlogis(1)=Inf
      initial <- bound(initial, 0.005)
      weights <- tmle_task$get_regression_task(update_node)$weights
      n <- length(unique(tmle_task$id))
      if(self$one_dimensional & for_fitting){
        # This computes (possibly weighted) ED and handles long case
        ED <- colSums(IC * weights)/n #apply(IC , 2, function(v) {sum(as.vector(matrix(v, nrow = n, byrow = T)*weights))})/length(weights)
      } else {
        ED <- NULL
      }

      if(length(observed) != length(initial)) {
        ratio <- length(initial) / length(observed)
        if(ratio%%1 == 0){
          warning("Observed and initial length do not match but are multiples of each other. Recycling values...")
          observed <- rep(observed, ratio)
        }
      }

      if(length(weights) != length(initial)) {
        ratio <- length(initial) / length(weights)
        if(ratio%%1 == 0){
          # This is for likelihood factors that output long_format predictions that dont match nrow of input task
          warning("Weights and initial length do not match but are multiples of each other. Recycling values...")
          weights <- rep(weights, ratio)
        }
      }

      submodel_data <- list(
        observed = observed,
        H = covariates_dt,
        initial = initial,
        submodel_info = submodel_info,
        ED = ED,
        update_node = update_node,
        weights = weights
      )



      if (drop_censored) {
        censoring_node <- tmle_task$npsem[[update_node]]$censoring_node$name
        if (!is.null(censoring_node)) {
          observed_node <- tmle_task$get_tmle_node(censoring_node)
          subset <- which(observed_node == 1)
          subset <- intersect(subset, which(tmle_task$get_tmle_node(update_node, compute_risk_set = T)[, at_risk] == 1))
          submodel_data <- list(
            observed = submodel_data$observed[subset],
            H = submodel_data$H[subset, , drop = FALSE],
            initial = submodel_data$initial[subset],
            submodel_info = submodel_info,
            ED = ED,
            update_node = update_node,
            weights = submodel_data$weights[subset]
          )
        }
      }

      return(submodel_data)
    },
    fit_submodel = function(submodel_data) {
      update_node <- submodel_data$update_node
      submodel_data["update_node"] <- NULL
      weights <- submodel_data$weights
      # TODO

      if(self$one_dimensional){
        # Will break if not called by original training task

        if(is.null(submodel_data$ED)) {
          #warning("No ED given in clever covariates. Defaulting to full EIC ED, which is incorrect.")
          submodel_data$H <- self$collapse_covariates(self$current_estimates, submodel_data$H)
          ED <- ED_from_estimates(self$current_estimates)
          EDnormed <- ED / (norm(ED, type = "2") )
          ED <- EDnormed


        } else {
          ED <- submodel_data$ED
          initial_variances <- private$.initial_variances

          vars <- unlist(lapply(initial_variances, `[[`, update_node))
          if(self$convergence_type == "scaled_var" & !is.null(vars)){
            #max_var <- max(vars)
            #Ensure that params with very small variances dont get too much weight
            median_var <- median(vars)
            min_var_allowed <- max(median_var/100,1e-3)
            zero <- vars < min_var_allowed
            vars[zero] <- min_var_allowed
            ED <- ED / sqrt(vars)
          }

          EDnormed <- ED / norm(ED, type = "2")# / sqrt(length(ED))))
          submodel_data$H <- submodel_data$H %*% EDnormed

          ED <- EDnormed


        }
      }

      submodel_data["ED"] <- NULL
      submodel_info <- submodel_data$submodel_info
      sub_index <- which(names(submodel_data) == "submodel_info")

      if (self$constrain_step) {
        ncol_H <- ncol(submodel_data$H)
        if (!(is.null(ncol_H) || (ncol_H == 1))) {
          stop(
            "Updater detected `constrain_step=TRUE` but multi-epsilon submodel.\n",
            "Consider setting `one_dimensional=TRUE`"
          )
        }


        risk <- function(epsilon) {

          submodel_estimate <- self$apply_submodel(submodel_data, epsilon)

          loss_function <- submodel_info$loss_function

          loss <- loss_function(submodel_estimate, submodel_data$observed) * weights
          mean(loss)

        }


        if (self$optim_delta_epsilon) {
          delta_epsilon <- self$delta_epsilon
          if(is.list(delta_epsilon)) {
            delta_epsilon <- delta_epsilon[[update_node]]
          }
          if(is.function(delta_epsilon)) {
            delta_epsilon <- delta_epsilon(submodel_data$H)
          }
          delta_epsilon <- c(0,delta_epsilon)

          min_eps = min(delta_epsilon)
          max_eps = max(delta_epsilon)

          optim_fit <- optim(
            par = list(epsilon = max_eps), fn = risk,
            lower = min_eps, upper = max_eps,
            method = "Brent"
          )
          epsilon <- optim_fit$par

        } else {
          epsilon <- self$delta_epsilon
          if(is.list(epsilon)) {
            epsilon <- epsilon[[update_node]]
          }
        }

        risk_val <- risk(epsilon)
        risk_zero <- risk(0)

         #TODO: consider if we should do this
        if(risk_zero<=risk_val){

          epsilon <- 0
          #private$.delta_epsilon <- private$.delta_epsilon/2
        }

        if (self$verbose) {
          cat(sprintf("risk_change: %e ", risk_val - risk_zero))
        }
      } else {
        if (self$fluctuation_type == "standard") {
          suppressWarnings({
            submodel_fit <- glm(observed ~ H - 1, submodel_data[-sub_index],
                                offset = submodel_info$offset_tranform(submodel_data$initial),
                                family = submodel_info$family,
                                weights = weights,

                                start = rep(0, ncol(submodel_data$H))
            )
          })
        } else if (self$fluctuation_type == "weighted") {
          if (self$one_dimensional) {
            suppressWarnings({
              submodel_fit <- glm(observed ~ -1, submodel_data[-sub_index],
                                  offset =  submodel_info$offset_tranform(submodel_data$initial),
                                  family = submodel_info$family,
                                  weights = as.numeric(H) * weights,
                                  start = rep(0, ncol(submodel_data$H))
              )
            })
          } else {
            warning(
              "Updater detected `fluctuation_type='weighted'` but multi-epsilon submodel.\n",
              "This is incompatible. Proceeding with `fluctuation_type='standard'`."
            )
            suppressWarnings({
              submodel_fit <- glm(observed ~ H - 1, submodel_data[-sub_index],
                                  offset =  submodel_info$offset_tranform(submodel_data$initial),
                                  family = submodel_info$family,
                                  weights = weights,
                                  start = rep(0, ncol(submodel_data$H))
              )
            })
          }
        }
        epsilon <- coef(submodel_fit)

        # NOTE: this protects against collinear covariates
        # (which we don't care about, we just want an update)
        epsilon[is.na(epsilon)] <- 0
      }

      if (self$verbose) {
        max_eps <- epsilon[which.max(abs(epsilon))]
        cat(sprintf("(max) epsilon: %e ", max_eps))
      }

      if(self$one_dimensional){

        epsilon <- epsilon * ED
      }
      return(epsilon)
    },
    submodel = function(epsilon, initial, H) {
      plogis(qlogis(initial) + H %*% epsilon)
    },
    loss_function = function(estimate, observed) {
      -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
    },
    apply_submodel = function(submodel_data, epsilon) {
      submodel_data$submodel_info$submodel(epsilon, submodel_data$initial, submodel_data$H)
    },
    apply_update = function(tmle_task, likelihood, fold_number, new_epsilon, update_node) {

      # get submodel data for all nodes
      submodel_data <- self$generate_submodel_data(
        likelihood, tmle_task,
        fold_number, update_node,
        drop_censored = FALSE
      )

      updated_likelihood <- self$apply_submodel(submodel_data, new_epsilon)

      if (any(!is.finite(updated_likelihood))) {
        stop("Likelihood was updated to contain non-finite values.\n
             This is likely a result of unbounded likelihood factors")
      }
      # un-scale to handle bounded continuous
      updated_likelihood <- tmle_task$unscale(
        updated_likelihood,
        update_node
      )

      return(updated_likelihood)
    },
    check_convergence = function(tmle_task, fold_number = "full") {
      estimates <- self$current_estimates

      n <- length(unique(tmle_task$id))
      if (self$convergence_type == "scaled_var") {
        # NOTE: the point of this criterion is to avoid targeting in an overly
        #       aggressive manner, as we simply need check that the following
        #       condition is met |P_n D*| / SE(D*) =< max(1/log(n), 1/10)
        # TODO Use variance from first iteration as weights so criterion does not change
        IC <- do.call(cbind, lapply(estimates, `[[`, "IC"))
        # TODO colVars is wrong when using long format
        # TODO The below is a correction that should be correct for survival (assuming long format is stacked by vectors of time and not by person)
        se_Dstar <- sqrt(apply(IC, 2, function(v){
          # If long then make it a matrix
          if(length(v)!=n){
            v <- matrix(v, nrow = n, byrow = T)
            #Collapse EIC for each person by summing across time (this is correct for survival)
            v <- rowSums(v)
          }

          return(var(v))
        })/n)
        # Handle case where variance is 0 or very small for whatever reason
        ED_threshold <- pmax(se_Dstar / min(log(n), 10), 1/n)
      } else if (self$convergence_type == "sample_size") {
        ED_threshold <- 1 / n
      }

      # get |P_n D*| of any number of parameter estimates
      ED <- ED_from_estimates(estimates)
      # zero out any that are from nontargeted parameter components
      ED <- ED * private$.targeted_components
      current_step <- self$step_number

      private$.EDs[[current_step]] <- ED


      ED_criterion <- abs(ED)

      if (self$verbose) {
        cat(sprintf("max(abs(ED)): %e\n", max(ED_criterion)))
      }


      return(all(ED_criterion <= ED_threshold))
    },
    update_best = function(likelihood) {
      current_step <- self$step_number
      ED <- private$.EDs[[current_step]]
      ED_2_norm <- sqrt(sum(ED^2))
      if (ED_2_norm < private$.best_ED) {
        likelihood$cache$update_best()
        private$.best_ED <- ED_2_norm
      }
    },
    update = function(likelihood, tmle_task) {
      update_fold <- self$update_fold
      maxit <- private$.maxit

      # seed current estimates
      private$.current_estimates <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$estimates(tmle_task, update_fold)
      })

      if(FALSE) {
        clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
          tmle_param$clever_covariates(tmle_task, update_fold)})
        IC <- lapply(clever_covariates, `[[`, "IC")
        if(!is.null(IC[[1]])){
          n <- length(unique(tmle_task$id))
          IC_vars <- lapply(IC, function(IC) {
            out <- lapply(self$update_nodes, function(node) {
              weights <- tmle_task$get_regression_task(node)$weights
              apply(IC[[node]] * weights,2, function(v) {var(rowSums(matrix(v, nrow = n, byrow = T)))})
            } )
            names(out) <- self$update_nodes
            return(out)
          })
          private$.initial_variances <- IC_vars


        } else {
          n <- length(unique(tmle_task$id))
          IC <- lapply(private$.current_estimates, `[[`, "IC")
          IC_vars <- lapply(IC, function(IC) {
            weights <- tmle_task$get_regression_task(node)$weights
            IC_var <- apply(IC[[node]] * weights,2, function(v) {var(rowSums(matrix(v, nrow = n, byrow = T)))})
            IC_var <- lapply(self$update_nodes, function(node) {IC_var})
            names(IC_var) <- self$update_nodes
            return(IC_var)
          })
          private$.initial_variances <- IC_vars
        }
      }


      #private$.initial_variances <- lapply(private$.current_estimates, `[[`, "var_comps")

      for (steps in seq_len(maxit)) {
        self$update_step(likelihood, tmle_task, update_fold)

        # update estimates based on updated likelihood
        private$.current_estimates <- lapply(self$tmle_params, function(tmle_param) {
          tmle_param$estimates(tmle_task, update_fold)
        })

        if (self$check_convergence(tmle_task, update_fold)) {
          break
        }

        if (self$use_best) {
          self$update_best(likelihood)
        }
      }

      if (self$use_best) {
        self$update_best(likelihood)
        likelihood$cache$set_best()
      }
    },
    register_param = function(new_params) {
      if (inherits(new_params, "Param_base")) {
        new_params <- list(new_params)
      }
      private$.tmle_params <- c(private$.tmle_params, new_params)
      private$.targeted_components <- unlist(lapply(private$.tmle_params, `[[`, "targeted"))
      new_update_nodes <- unlist(lapply(new_params, `[[`, "update_nodes"))
      private$.update_nodes <- unique(c(
        private$.update_nodes,
        new_update_nodes
      ))
    },
    set_estimates = function(tmle_task, update_fold = "full"){
      private$.current_estimates <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$estimates(tmle_task, update_fold)
      })
      #TODO Variance weights should be based on each component separately.
      # Might be better to have estimates return IC and IC-components.
      # private$.initial_variances <- lapply(private$.current_estimates, function(ests) {
      #   resample::colVars(matrix(ests$IC, nrow = length(unique(tmle_task$id))))
      # })
      private$.initial_variances <- lapply(private$.current_estimates, `[[`, "var_comps")
    }
  ),
  active = list(
    epsilons = function() {
      return(private$.epsilons)
    },
    EDs = function() {
      return(private$.EDs)
    },
    tmle_params = function(new_params = NULL) {
      if (!is.null(new_params)) {
        if (inherits(new_params, "Param_base")) {
          new_params <- list(new_params)
        }
        private$.tmle_params <- new_params
        private$.update_nodes <- unique(unlist(lapply(
          new_params, `[[`,
          "update_nodes"
        )))
      }
      return(private$.tmle_params)
    },
    update_nodes = function() {
      return(private$.update_nodes)
    },
    update_fold = function() {
      if (self$cvtmle) {
        # use training predictions on validation sets
        update_fold <- "validation"
      } else {
        # use predictions from full fit
        update_fold <- "full"
      }
    },
    step_number = function() {
      return(private$.step_number)
    },
    maxit = function() {
      return(private$.maxit)
    },
    cvtmle = function() {
      return(private$.cvtmle)
    },
    one_dimensional = function() {
      return(private$.one_dimensional)
    },
    constrain_step = function() {
      return(private$.constrain_step)
    },
    delta_epsilon = function() {
      return(private$.delta_epsilon)
    },
    convergence_type = function() {
      return(private$.convergence_type)
    },
    fluctuation_type = function() {
      return(private$.fluctuation_type)
    },
    optim_delta_epsilon = function() {
      return(private$.optim_delta_epsilon)
    },
    use_best = function() {
      return(private$.use_best)
    },
    verbose = function() {
      return(private$.verbose)
    },
    current_estimates = function() {
      return(private$.current_estimates)
    },
    initial_variances = function(){
      return(private$.initial_variances)
    }
  ),
  private = list(
    .epsilons = list(),
    .EDs = list(),
    .best_ED = Inf,
    .tmle_params = NULL,
    .update_nodes = NULL,
    .step_number = 0,
    # TODO: change maxit for test
    .maxit = 100,
    .cvtmle = NULL,
    .one_dimensional = NULL,
    .constrain_step = NULL,
    .delta_epsilon = NULL,
    .optim_delta_epsilon = NULL,
    .convergence_type = NULL,
    .fluctuation_type = NULL,
    .use_best = NULL,
    .verbose = FALSE,
    .targeted_components = NULL,
    .current_estimates = NULL,
    .initial_variances = NULL
  )
)
