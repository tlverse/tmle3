#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_CR <- R6Class(
  classname = "Param_CR",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, censoring_node = "processA", competing_risk_nodes, target_risk_node, target_times = NULL) {
      super$initialize(observed_likelihood, list(), target_risk_node)
      private$.update_nodes <- union(competing_risk_nodes, target_risk_node)
      private$.cf_likelihood <- CF_Likelihood$new(observed_likelihood, intervention_list)
      private$.cf_tasks <- private$.cf_likelihood$cf_tasks
      private$.target_times <- target_times
      private.nodes <- list(competing_risk_nodes = setdiff(competing_risk_nodes, target_risk_node),
                            censoring_node = censoring_node,
                            target_risk_node = target_risk_node )

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # Get observed g likelihoods
      g <- self$likelihood$get_likelihoods(tmle_task, names(self$intervention_list), fold_number = fold_number)
      g_cf <- self$cf_likelihood$get_likelihoods(tmle_task, names(self$intervention_list), fold_number = fold_number)
      target_node <- self$competing_nodes$target_risk_node
      competing_risk_nodes <- self$competing_nodes$competing_risk_nodes
      censoring_node <-  self$competing_nodes$censoring_node
      all_nodes <- c(censoring_node, competing_risk_nodes, target_node)


      # assumes all competing risks have same times
      times <- unique(sort(tmle_task$npsem[[target_node]]$time))
      num_times <- length(times)

      # Those are still at risk
      at_risk_indicator <- tmle_task$get_tmle_node(target_node, expand = T, compute_risk_set = T)$at_risk
      at_risk_indicator_mat <- self$long_to_wide(at_risk_indicator_mat, num_times)

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
      #Ensures we get actual survival/hazard probabilities for each time (conditional that they are at risk) and not the degenerate values
      cf_task$force_at_risk <- T
      G <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$censoring_node, fold_number = fold_number)
      Gmat <- self$long_to_wide(G, num_times)
      Gmat <- cbind(rep(1, nrow(G)), Gmat[,-ncol(Gmat)])
      # Get competing risk hazards
      Qcompeting <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$competing_risk_nodes, fold_number = fold_number)
      Qtarget <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$target_risk_node, fold_number = fold_number)
      Qtarget_mat <- self$long_to_wide(Qtarget, num_times)
      observed <- lapply(c(target_node, ))
      #get cumulative hazard
      hazard_dt <- as.data.table(cbind(Qcompeting, Qtarget))
      setnames(hazard_dt, c(self$competing_nodes$competing_risk_node,  self$competing_nodes$target_risk_node))
      cum_hazard <- self$cumulative_hazard(hazard_dt)
      cum_hazard_mat <- self$long_to_wide(hazard_dt, num_times)
      cum_survival <- self$hazard_to_survival(cum_hazard)

      cum_inc_mat <- self$conditional_cumulative_indicence(cum_survival, Qtarget_mat)
      target_times <- self$target_times
      # Only theose still with nondegenerate values are updated
      Hg = g_cf/(g*Gmat)*at_risk_indicator_mat
      #long format but filled byrow to match predictions
      H_target <- do.call(cbind, lapply(target_times, function(t_tgt) {
        res <- Hg*(1-(cum_inc_mat[,t_tgt] - cum_inc_mat))/cum_survival
        res[, seq_along(res) > t_tgt] <- 0
        return(as.vector(t(res)))
      }))
      H_competitors <- do.call(cbind, lapply(target_times, function(t_tgt) {
        res <- -1 * Hg*(cum_inc_mat[,t_tgt] - cum_inc_mat)/cum_survival
        res[, seq_along(res) >= t_tgt] <- 0
        return(as.vector(t(res)))
      }))
      H_list <- list()
      H_list[[target_node]] <- H_target
      for(node in competing_risk_nodes) {
        H_list[[node]] <- H_competitors
      }

      if(tmle_task$uuid == self$observed_likelihood$uuid){
        #If training task then compute and add the EIC mean for one-step
        D_list <- list()

        obs_vals <- unlist(tmle_task$get_tmle_node(target_node, format = T, expand = T, compute_risk_set = F), use.names = F)
        Qvals <- Qtarget
        residuals <- (obs_vals - Qvals)
        D_list[[node]] <- colSums(H_target*residuals)/nrow(cum_hazard_mat)

        D_list[[target_node]] <- D_tgt
        for(node in competing_risk_nodes) {
          obs_vals <- unlist(tmle_task$get_tmle_node(node, format = T, expand = T, compute_risk_set = F), use.names = F)
          Qvals <- Qcompeting[[node]]
          residuals <- (obs_vals - Qvals)
          D_list[[node]] <- colSums(H_competitors*residuals)/nrow(cum_hazard_mat)
        }
        H_list[["ED"]] <- D_list
      }

      return(H_list)
    },
    cumulative_hazard = function(hazard_dt){
      collapsed_hazard <- colSums(hazard_dt)
      return(as.vector(collapsed_hazard))
    },
    long_to_wide = function(long_vec, num_times, byrow = T){
      matrix(long_vec, ncol = num_times, byrow = byrow)
    },
    hazard_to_survival = function(hazard_mat){
      t(apply(1-hazard_mat,1,cumprod))
    },
    conditional_cumulative_indicence_mat = function(cum_survival_mat, hazard_target_mat){
      result <- matrix(0, nrow = nrow(cum_survival_mat), ncol = ncol(cum_survival_mat))
      for(t in 1:ncol(result)) {
        if(t==1){
          result[,t] <- hazard_target_mat[,1]
        } else {
          result <- result[,t - 1] + cum_survival_mat[, t-1]*hazard_target_mat[,t]
        }
      }
      return(result)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      # Get observed g likelihoods
      g <- self$likelihood$get_likelihoods(tmle_task, names(self$intervention_list), fold_number = fold_number)
      g_cf <- self$cf_likelihood$get_likelihoods(tmle_task, names(self$intervention_list), fold_number = fold_number)
      target_node <- self$competing_nodes$target_risk_node
      competing_risk_nodes <- self$competing_nodes$competing_risk_nodes
      censoring_node <-  self$competing_nodes$censoring_node
      all_nodes <- c(censoring_node, competing_risk_nodes, target_node)


      # assumes all competing risks have same times
      times <- unique(sort(tmle_task$npsem[[target_node]]$time))
      num_times <- length(times)

      # Those are still at risk
      at_risk_indicator <- tmle_task$get_tmle_node(target_node, expand = T, compute_risk_set = T)$at_risk
      at_risk_indicator_mat <- self$long_to_wide(at_risk_indicator_mat, num_times)

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
      #Ensures we get actual survival/hazard probabilities for each time (conditional that they are at risk) and not the degenerate values
      cf_task$force_at_risk <- T
      G <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$censoring_node, fold_number = fold_number)
      Gmat <- self$long_to_wide(G, num_times)
      Gmat <- cbind(rep(1, nrow(G)), Gmat[,-ncol(Gmat)])
      # Get competing risk hazards
      Qcompeting <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$competing_risk_nodes, fold_number = fold_number)
      Qtarget <- self$likelihood$get_likelihoods(cf_task, self$competing_nodes$target_risk_node, fold_number = fold_number)
      Qtarget_mat <- self$long_to_wide(Qtarget, num_times)
      observed <- lapply(c(target_node, ))
      #get cumulative hazard
      hazard_dt <- as.data.table(cbind(Qcompeting, Qtarget))
      setnames(hazard_dt, c(self$competing_nodes$competing_risk_node,  self$competing_nodes$target_risk_node))
      cum_hazard <- self$cumulative_hazard(hazard_dt)
      cum_hazard_mat <- self$long_to_wide(hazard_dt, num_times)
      cum_survival <- self$hazard_to_survival(cum_hazard)

      cum_inc_mat <- self$conditional_cumulative_indicence(cum_survival, Qtarget_mat)
      target_times <- self$target_times
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },

    update_nodes = function() {
      return(c(private$.update_nodes))
    },
    gradient = function(){
      private$.gradient
    },
    competing_nodes = function(){
      private$.nodes
    },
    target_times = function(){
      private$.target_times
    }
  ),
  private = list(
    .type = "ATE",
    .cf_likelihood = NULL,
    .cf_tasks = NULL,
    .supports_outcome_censoring = FALSE,
    .submodel_type_supported = c("logistic"),
    .update_nodes = NULL,
    .nodes = NULL,
    .target_times = NULL
  )
)
