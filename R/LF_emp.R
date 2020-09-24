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

      # TODO weights
      observed <- (tmle_task$get_tmle_node(self$name, format = T, include_id = F))
      if(ncol(observed) == 0 | length(observed) ==0){
        private$.empirical_fit <- list(emp_probs = NULL, uniq_obs = NULL )
        return()

      }
      uniq_obs <- unique(observed)

      match_index <- uniq_obs[observed, which = T, on = colnames(uniq_obs)]
      weights_mat <- data.table(weights = weights, grp = match_index)
      emp_probs <- weights_mat[, sum(weights), by = grp]
      emp_probs$grp <- NULL
      emp_probs <- unlist(emp_probs) / sum(weights)
      if(length(emp_probs) != nrow(uniq_obs)) {
        stop("LF_emp error")
      }
      #emp_probs <- table(match_index) /nrow(observed)

      #SUPER SLOW. Thank you data.table
      # counts <- unlist(apply(uniq_obs, 1, function(obs){ sum(weights * as.numeric(apply(observed, 1, function(new_obs){
      #
      #   all(new_obs == obs)})))}))
      #
      # emp_probs <- counts/sum(counts)


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
      observedfull <-  tmle_task$get_tmle_node(self$name, format = T, include_id = T, include_time = T, expand = expand)
      if(length(observedfull)==0){
        return(observedfull)
      }

      observed <- observedfull[, node, with = F] #unlist(observedfull[, setdiff(colnames(observedfull), c("id", "t")), with = F])

      # TODO dont need weights for prediction??
      #weights <- tmle_task$get_regression_task(self$name)$weights
      emp_probs <-  private$.empirical_fit$emp_probs
      if(is.null(emp_probs)) {
        return(rep(1/tmle_task$nrow, tmle_task$nrow))
      }
      uniq_obs <-  private$.empirical_fit$uniq_obs

      match_index <- uniq_obs[observed, which = T, on = colnames(uniq_obs)]
      probs <- emp_probs[match_index]
      probs[is.na(probs)] <- 0

      # TODO THe above might be wrong.
      # matched_obs <- apply(observed, 1, function(obs){
      #   index <- which(unlist(apply(uniq_obs, 1, function(uniq) {
      #     all(uniq == obs)
      #   })))
      #   if(length(index) == 0){
      #     return(NA)
      #   }
      #   return(index)
      # })
       # match(observed, uniq_obs)
      #Match observations to empirical probs obtained from training set


      #probs <- weights*probs
      #weights <- tmle_task$get_regression_task(self$name)$weights
      #probs <- as.data.table(probs )
      #setnames(probs, self$name)
      #probs$id = observedfull$id
      #probs$t = observedfull$t
      probs <- unlist(probs)
      # probs$t = NULL

      return(probs)
    },
    sample = function(tmle_task = NULL, n_samples = NULL, fold_number = "full") {
      #TODO No conditioning variables so dont need task
      if(is.null(tmle_task)) {
        num <- 1
      } else {
        observedfull <-  tmle_task$get_tmle_node(self$name, format = T, include_id = T, include_time = T)
        num <- 1:nrow(observedfull)
      }
      emp_probs <-  private$.empirical_fit$emp_probs
      uniq_obs <-  private$.empirical_fit$uniq_obs

      values <- as.matrix(do.call(rbind, lapply(1:num, function(i) uniq_obs[base::sample(nrow(uniq_obs), n_samples, prob = emp_probs, replace = T)])))
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
