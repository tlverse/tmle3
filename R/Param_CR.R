#' Competing Risks cause-specific cumulative incidence (no ties)
#' Supports multiple intervention nodes, and multiple competing risks + censoring.
#' Parameter definition for the CR cumulative incidence.
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#' @export
Param_CR <- R6Class(
  classname = "Param_CR",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, marginalized = F, censoring_node = "processA", competing_risk_nodes, target_risk_node, target_times = NULL, node_time_ordering = NULL) {
      if(!marginalized & is.null(node_time_ordering)) {
        stop("You specified to target via the non-marginalized cause specific hazards but did provide the time ordering of the regression.")
      }
      private$.update_nodes <- union(competing_risk_nodes, target_risk_node)
      super$initialize(observed_likelihood, list(), outcome_node = target_risk_node)

      private$.cf_likelihood <- CF_Likelihood$new(observed_likelihood, intervention_list)
      private$.cf_tasks <- private$.cf_likelihood$cf_tasks
      private$.target_times <- target_times
      private$.marginalized <- marginalized
      private$.nodes <- list(competing_risk_nodes = setdiff(competing_risk_nodes, target_risk_node),
                            censoring_node = censoring_node,
                            target_risk_node = target_risk_node, node_time_ordering = node_time_ordering )

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_estimates = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # Get observed g likelihoods
      intervention_nodes <- names(self$intervention_list)

      g <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number = fold_number)
      g_cf <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number = fold_number)

      # collapse across multiple intervention nodes
      if (!is.null(ncol(g)) && ncol(g) > 1) {
        g <- apply(g,1, prod)
        g <- bound(g, c(0.005,1))
      }
      # collapse across multiple intervention nodes
      if (!is.null(ncol(g_cf)) && ncol(g_cf) > 1) {
        g_cf <- apply(g_cf,1, prod)
      }

      target_node <- self$competing_nodes$target_risk_node
      competing_risk_nodes <- self$competing_nodes$competing_risk_nodes
      censoring_node <-  self$competing_nodes$censoring_node
      all_nodes <- c(censoring_node, competing_risk_nodes, target_node)


      # assumes all competing risks have same times
      times <- unique(sort(tmle_task$npsem[[target_node]]$time))
      num_times <- length(times)

      # Those are still at risk
      at_risk_indicator <- tmle_task$get_tmle_node(target_node, expand = T, compute_risk_set = T)$at_risk
      at_risk_indicator_mat <- self$long_to_wide(at_risk_indicator, num_times)

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
      #Ensures we get actual survival/hazard probabilities for each time (conditional that they are at risk) and not the degenerate values
      cf_task$force_at_risk <- T
      G <- self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$censoring_node, fold_number = fold_number)
      G <- bound(G, c(0.005,1))
      Gmat <- self$long_to_wide(G, num_times)
      Gmat <- self$hazard_to_survival(Gmat)
      Gmat <- cbind(rep(1, nrow(Gmat)), Gmat[,-ncol(Gmat)])
      # Get competing risk hazards
      Qcompeting <- as.matrix(self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$competing_risk_nodes, fold_number = fold_number))
      Qcompeting <- bound(Qcompeting, c(0.005,1))
      colnames(Qcompeting) <- competing_risk_nodes
      Qtarget <- unlist(self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$target_risk_node, fold_number = fold_number))
      Qtarget <- bound(Qtarget, c(0.005,1))
      Qtarget_mat <- self$long_to_wide((Qtarget), num_times)

      #get cumulative hazard
      hazard_dt <- as.matrix(cbind(Qcompeting, Qtarget))
      colnames(hazard_dt) <-  c(self$competing_nodes$competing_risk_node,  self$competing_nodes$target_risk_node)

      if(!self$marginalized) {
        #If targetting is not through marginalized hazards then compute marginalized hazard
        ordering <- unlist(self$competing_nodes$node_time_ordering)
        hazard_dt <- hazard_dt[, ordering]
        expand_dt <- t(apply(1 - hazard_dt, 1 , cumprod))

        expand_dt <- cbind(rep(1, nrow(expand_dt)), expand_dt[, -ncol(expand_dt)])

        hazard_dt <- expand_dt * hazard_dt
        colnames(hazard_dt) <- ordering
        marginal_Qtgt_mat <- self$long_to_wide(hazard_dt[, target_node], num_times)

      } else {
        marginal_Qtgt_mat <- Qtarget_mat
      }

      #cumulative is hazard is the sum of the marginalized hazards
      cum_hazard <- self$cumulative_hazard(hazard_dt)

      cum_hazard_mat <- self$long_to_wide(cum_hazard, num_times)

      cum_survival <- self$hazard_to_survival(cum_hazard_mat)
      cum_inc_mat <- self$conditional_cumulative_indicence_mat(cum_survival, marginal_Qtgt_mat)
      target_times <- self$target_times
      if(is.null(target_times)){
        target_times <- times
        private$.target_times <- target_times
      }
      # Only theose still with nondegenerate values are updated
      Hg = g_cf/(g*Gmat) #* (at_risk_indicator_mat)
      #long format but filled byrow to match predictions
      H_target <- do.call(cbind, lapply(target_times, function(t_tgt) {
        res <- Hg*(1-(cum_inc_mat[,t_tgt] - cum_inc_mat))/cum_survival
        res[, seq_len(ncol(res)) > t_tgt] <- 0
        return(bound(as.vector(t(res)), c(-40,40)))
      }))
      H_competitors <- do.call(cbind, lapply(target_times, function(t_tgt) {
        res <- -1 * Hg*(cum_inc_mat[,t_tgt] - cum_inc_mat)/cum_survival
        res[, seq_len(ncol(res)) >= t_tgt] <- 0
        return(bound(as.vector(t(res)), c(-40,40)))
      }))
      H_list <- list()
      H_list[[target_node]] <- H_target *as.vector(at_risk_indicator)

      for(node in competing_risk_nodes) {
        at_risk_indicator <- tmle_task$get_tmle_node(node, expand = T, compute_risk_set = T)$at_risk

        H_list[[node]] <- H_competitors * at_risk_indicator

      }

      if(tmle_task$uuid == self$observed_likelihood$training_task$uuid){
        #If training task then compute and add the EIC mean for one-step
        D_list <- list()
        EIC_list <- list()
        obs_vals <- unlist(tmle_task$get_tmle_node(target_node, format = T, expand = T, compute_risk_set = F), use.names = F)
        Qvals <- Qtarget
        residuals <- as.vector(obs_vals - Qvals)

        EIC_tgt <- H_list[[target_node]]*residuals*tmle_task$get_regression_task(target_node)$weights

        if(for_estimates) {
          EIC_list[[target_node]] <- EIC_tgt
        }

        D_list[[target_node]] <- colSums(EIC_tgt)/nrow(cum_hazard_mat)

        for(node in competing_risk_nodes) {
          obs_vals <- unlist(tmle_task$get_tmle_node(node, format = T, expand = T, compute_risk_set = F), use.names = F)
          Qvals <- Qcompeting[,node]
          residuals <- as.vector(obs_vals - Qvals)
          EIC_comp <- H_list[[node]]*residuals*tmle_task$get_regression_task(node)$weights
          if(for_estimates) {
            EIC_list[[node]] <- EIC_comp
          }

          D_list[[node]] <- colSums(EIC_comp)/nrow(cum_hazard_mat)
        }
        H_list[["ED"]] <- D_list

        if(for_estimates) {
          H_list[["EIC"]] <- EIC_list
        }

      }

      return(H_list)
    },
    cumulative_hazard = function(hazard_dt){
      collapsed_hazard <- as.vector(apply(hazard_dt, 1, sum))
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
          result[,t] <- result[,t - 1] + cum_survival_mat[, t-1]*hazard_target_mat[,t]
        }
      }
      return(result)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      # Get observed g likelihoods
      target_node <- self$competing_nodes$target_risk_node
      competing_risk_nodes <- self$competing_nodes$competing_risk_nodes
      censoring_node <-  self$competing_nodes$censoring_node
      all_nodes <- c(censoring_node, competing_risk_nodes, target_node)


      # assumes all competing risks have same times
      times <- unique(sort(tmle_task$npsem[[target_node]]$time))
      num_times <- length(times)

      # Those are still at risk
      at_risk_indicator <- tmle_task$get_tmle_node(target_node, expand = T, compute_risk_set = T)$at_risk
      at_risk_indicator_mat <- self$long_to_wide(at_risk_indicator, num_times)

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
      #Ensures we get actual survival/hazard probabilities for each time (conditional that they are at risk) and not the degenerate values
      cf_task$force_at_risk <- T
      G <- self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$censoring_node, fold_number = fold_number)
      G <- bound(G, c(0.005,1))
      Gmat <- self$long_to_wide(G, num_times)
      Gmat <- self$hazard_to_survival(Gmat)
      Gmat <- cbind(rep(1, nrow(Gmat)), Gmat[,-ncol(Gmat)])
      # Get competing risk hazards
      Qcompeting <- self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$competing_risk_nodes, fold_number = fold_number)
      Qcompeting <- bound(Qcompeting, c(0.005,1))
      Qtarget <- self$observed_likelihood$get_likelihoods(cf_task, self$competing_nodes$target_risk_node, fold_number = fold_number)
      Qtarget <- bound(Qtarget, c(0.005,1))

      Qtarget_mat <- self$long_to_wide(Qtarget, num_times)

      #get cumulative hazard
      hazard_dt <- as.matrix(cbind(Qcompeting, Qtarget))
      colnames(hazard_dt) <-  c(self$competing_nodes$competing_risk_node,  self$competing_nodes$target_risk_node)

      if(!self$marginalized) {
        #If targetting is not through marginalized hazards then compute marginalized hazard
        ordering <- unlist(self$competing_nodes$node_time_ordering)
        hazard_dt <- hazard_dt[, ordering]
        expand_dt <- t(apply(1 - hazard_dt, 1 , cumprod))

        expand_dt <- cbind(rep(1, nrow(expand_dt)), expand_dt[, -ncol(expand_dt)])

        hazard_dt <- expand_dt * hazard_dt
        colnames(hazard_dt) <- ordering
        marginal_Qtgt_mat <- self$long_to_wide(hazard_dt[, target_node], num_times)

      } else {
        marginal_Qtgt_mat <- Qtarget_mat
      }

      cum_hazard <- self$cumulative_hazard(hazard_dt)

      cum_hazard_mat <- self$long_to_wide(cum_hazard, num_times)

      cum_survival <- self$hazard_to_survival(cum_hazard_mat)
      cum_inc_mat <- self$conditional_cumulative_indicence_mat(cum_survival, marginal_Qtgt_mat)

      clevs <- self$clever_covariates(tmle_task, fold_number, for_estimates = T)
      target_times <- self$target_times

      EIC <- lapply(clevs$EIC, function(X) {
          apply(X,2, function(v) {
        rowSums(matrix(v, nrow = nrow(cum_inc_mat), byrow = T))
      } ) })
      var_comps <- lapply(EIC, resample::colVars)
      names(var_comps) <- names(clevs$EIC)
      EIC <- Reduce('+', EIC) + t(t(cum_inc_mat[, target_times]) - as.vector(colMeans(cum_inc_mat[, target_times])))
      return(list(IC = EIC, var_comps = var_comps, psi = colMeans(cum_inc_mat[, target_times])))
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E_W P[T<=t, K = %s| %s, W]", self$outcome_node, paste(self$cf_likelihood$name, collapse = ", ") )
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
    },
    marginalized = function(){
      return(private$.marginalized)
    }
  ),
  private = list(
    .type = "CR",
    .cf_likelihood = NULL,
    .cf_tasks = NULL,
    .supports_outcome_censoring = FALSE,
    .submodel_type_supported = c("logistic"),
    .update_nodes = NULL,
    .nodes = NULL,
    .target_times = NULL,
    .marginalized = NULL
  )
)
