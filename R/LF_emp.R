#' Likelihood Factor Estimated using Empirical Distribution
#'
#' Uses the empirical probability distribution (puts mass \eqn{1/n} on each of the observations, or uses weights if specified) to estimate a marginal density.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#' Only compatible with marginal likelihoods (no parent nodes). Only compatible with densities (no conditional means).
#' The \code{type} argument will be ignored if specified.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_emp, name, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @export
LF_emp <- R6Class(
  classname = "Lf_emp",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, ...) {
      super$initialize(name, ..., type = "density")
      private$.name <- name
    },
    train = function(tmle_task, ...) {
      # get possible values from task if discrete
      super$train(tmle_task)
      weights <- tmle_task$get_regression_task(self$name)$weights

      observed <- unlist(tmle_task$get_tmle_node(self$name, format = T, include_id = F), use.names = F)
      uniq_obs <- unique(observed)
      counts <- unlist(lapply(uniq_obs, function(obs){ sum(weights * as.numeric(observed == obs))}), use.names = F)

      emp_probs <- counts/sum(counts)
      names(emp_probs) <- uniq_obs

      private$.empirical_fit <- list(emp_probs = emp_probs, uniq_obs = uniq_obs )
    },
    get_mean = function(tmle_task, fold_number = "full") {
      stop("nothing to predict")
    },
    get_density = function(tmle_task, fold_number = "full", expand = T, node = NULL) {
      if (is.null(node)) {
        node <- tmle_task$npsem[[self$name]]$variables
      } else {
        node <- tmle_task$npsem[[node]]$variables
      }
      # TODO: this only makes sense if the tmle_task is the same as the training one
      # Does not use the empirical measure of training sample but instead recomputes each time for new data
      # Not sure if this is what we want.
      # TODO: this only makes sense if the tmle_task is the same as the training one
      # This computes the true empirical density, including when there are ties.
      observedfull <-  tmle_task$get_tmle_node(self$name, format = T, include_id = T, include_time = T, expand = expand)


      observed <-observedfull[, node, with = F][[1]] #unlist(observedfull[, setdiff(colnames(observedfull), c("id", "t")), with = F])

      # TODO dont need weights for prediction??
      #weights <- tmle_task$get_regression_task(self$name)$weights
      emp_probs <-  private$.empirical_fit$emp_probs
      uniq_obs <-  private$.empirical_fit$uniq_obs


      matched_obs <- match(observed, uniq_obs)
      #Match observations to empirical probs obtained from training set
      probs <- as.vector(sapply(matched_obs, function(i) {
        if(is.na(i)) {
          return(0)
        }
        return(emp_probs[i])
      }))

      #probs <- weights*probs
      #weights <- tmle_task$get_regression_task(self$name)$weights
      probs <- data.table(probs )
      setnames(probs, self$name)
      probs$id = observedfull$id
      probs$t = observedfull$t

      # probs$t = NULL
      return(probs)
    },
    sample = function(tmle_task = NULL, n_samples = NULL, fold_number = "full") {
      #TODO No conditioning variables so dont need task
      observedfull <-  tmle_task$get_tmle_node(self$name, format = T, include_id = T, include_time = T)

      emp_probs <-  private$.empirical_fit$emp_probs
      uniq_obs <-  private$.empirical_fit$uniq_obs

      values <- as.matrix(do.call(rbind, lapply(1:nrow(observedfull), function(i) uniq_obs[base::sample(seq_along(uniq_obs), n_samples, prob = emp_probs)])))
      return(values)
    }
  ),
  active = list(
    empirical_fit = function()
      private$.empirical_fit
  ),
  private = list(
    .name = NULL,
    .empirical_fit = NULL
  )
)
