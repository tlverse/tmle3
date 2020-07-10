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
                              fit_method = "l2") {
      super$initialize(maxit = maxit, cvtmle=cvtmle, 
                       one_dimensional = one_dimensional,
                       constrain_step = constrain_step , 
                       delta_epsilon = delta_epsilon,
                       convergence_type = convergence_type,
                       fluctuation_type = fluctuation_type,
                       use_best = use_best,
                       verbose = verbose)
     private$.fit_method <- fit_method
    },
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
            ind <- max(which(norms <= self$delta_epsilon))
            if (ind > 1) break
            
            fit_lambda <- submodel_fit$lambda
            
            if(fit_lambda==1){
              stop("only one lambda could be fit")
            }
            
            # try to estimate what the correct lambda value is and go a bit beyond that
            norm_ratio <- self$delta_epsilon/norms[2]
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
          max_eps <- epsilon[which.max(abs(epsilon))]
          cat(sprintf("(max) epsilon: %e ", max_eps))
        }
      } else {
        # TODO: check
        # print("classic")
        epsilon <- super$fit_submodel(submodel_data)
      } 
    
      return(epsilon)
    }
  ),
  active = list(
    fit_method = function() {
      return(private$.fit_method)
    }

  ),
  private = list(
    .fit_method = NULL
  )
)
