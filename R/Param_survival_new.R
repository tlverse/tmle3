#' Survival Curve
#'
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
#'   \code{define_param(Param_survival, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_survival <- R6Class(
  classname = "Param_survival",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ..., outcome_node, target_times = NULL) {

       # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      # W A processN processA

      private$.cf_likelihood <- CF_Likelihood$new(observed_likelihood, intervention_list)
      times <- setdiff(sort(unique(observed_likelihood$training_task$time)),0)
      private$.times <- times
      if(is.null(target_times)){
        target_times <- times
        private$.targeted <- rep(TRUE, length(target_times))
      } else {
        private$.targeted <-  rep(TRUE, length(target_times))
      }
      private$.target_times <- target_times

      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
    },
    long_to_mat = function(x,id, time){
      dt <- data.table(id=id,time=time,x=as.vector(x))
      wide <- dcast(dt, id~time, value.var="x")
      mat <- as.matrix(wide[,-1,with=FALSE])
      return(mat)
    },
    hm_to_sm = function(hm){
      # TODO: check
      sm <- t(apply(1-hm,1,cumprod))
      # sm <- cbind(1,sm[,-ncol(sm)])
      return(sm)
    },
    clever_covariates_internal = function(tmle_task = NULL, fold_number = "full", subset_times = FALSE) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      is_observed <- tmle_task$uuid ==  self$observed_likelihood$training_task$uuid
      # Intervene on censoring/riskset so that hazard probabilities are returned
      # In particular, set the counting processes to 1 and specify that risk set should not be checked


      cf_task = task$clone()
      cf_task$force_at_risk <- TRUE
      #cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
      intervention_nodes <- names(self$intervention_list)

      g_A <- self$observed_likelihood$get_likelihoods(cf_task, "A", fold_number, drop_id = T)
      #g_A <- unlist(g_A[,"A", with = F])
      g_A <- bound(g_A, c(0.005,1))
      # I(A=1)
      #TODO this
      cf_pA <- self$cf_likelihood$get_likelihoods(cf_task, "A", fold_number, drop_id = T)

      #cf_pA <- unlist(cf_pA[,"A", with = F], use.names = F)

      #cf_pA <- tmle_task$get_tmle_node("A")[,A]

      #Should already be matrix. Assumed to be in order of time
      Q <- self$observed_likelihood$get_likelihoods(cf_task, c("processN"), fold_number)
      # TODO: make bound configurable
      n = length(cf_pA)

      Q <- matrix(Q, nrow = n, byrow = T)

      Q <- bound(Q, c(0.005,1))


      G <- self$observed_likelihood$get_likelihoods(cf_task, c("processA"), fold_number)

      G <- matrix(G, nrow = n, byrow = T)

      G <- bound(G, c(0.005,1))

      Q_surv <- as.matrix(self$hm_to_sm(Q))
      G_surv <- as.matrix(self$hm_to_sm(G))

      # fix t-1
      G_surv <- cbind(rep(1, nrow(G_surv)),G_surv[,-ncol(G_surv)])

      time <- cf_task$data$t

      times <- setdiff(sort(unique(time)),0)
      #t_mat = matrix(rep(ks,nrow(cf_pA)), nrow = nrow(cf_pA), ncol = length(ks), byrow =T)

      # Compute clever covariates for one  tgt time in long format
      # stacked by vertically in chunks of time

      hk_at_tgt <- function(tgt_time) {
        index_set_to_zero <- which(times > tgt_time)


        # Get survival at target time
        Q_surv_tgt <- Q_surv[,which(tgt_time ==times)]
        # This is the clever covariate matrix at a single target time


        clever_cov = -1 *(cf_pA*Q_surv_tgt/g_A)/ (G_surv*Q_surv)
        if(length(index_set_to_zero) >0){
          clever_cov[,index_set_to_zero] <- 0
        }

        return(as.vector(t(clever_cov)))
      }
      target_times <- self$target_times
      # bind the clever covariates for each tgt time by columns
      HA <- as.matrix(do.call(cbind, lapply(target_times, hk_at_tgt)))

      # get observed counting process
      # Contains columns t and id

      observed_N <- tmle_task$get_tmle_node("processN", include_time = T, include_id = T, expand= T)
      at_risk <- tmle_task$get_regression_task("processN", drop_censored = F, expand = T)$get_data(,"at_risk")

      # Only those at risk have likelihood updated.


      if(is_observed){
        #TODO check
        observed_N_wide <- reshape(observed_N, idvar = "id", timevar = "t", direction = "wide")
        observed_N_wide$id <- NULL
        # TODO this is crazy slow, code better and dont recompute
        to_dNt <- function(v){

          dt = c(0,diff((v)))
          return(dt)
        }
        # Converts to dNt format
        observed_dN_wide <- t(apply(observed_N_wide, 1, to_dNt))




        jump_time <- nrow(observed_N_wide) - apply(observed_N_wide, 1, sum) + 1
        t_mat <- matrix(1:ncol(observed_dN_wide), ncol = ncol(observed_dN_wide), nrow = nrow(observed_dN_wide), byrow = T )
        # TODO dont recompute this every time
        ind = t_mat <= jump_time



        # compute scaled and zeroed residuals vector


        residuals = as.vector(t((as.matrix(observed_dN_wide) - Q)*ind))

        clever_dot_HA <- HA*residuals


        private$.D_cache[[tmle_task$uuid]] <-  clever_dot_HA



      }
      zero_rows <- which(at_risk == 0)
      # TODO dont compute HA for all people
      HA[zero_rows,] <- 0

      return(list(processN =  HA))
    },
    get_EIC_var = function(tmle_task, fold_number = "full"){
      self$clever_covariates_internal(tmle_task, fold_number)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      #TODO compute estimates and EIC
      # computes counterfactual survival curve via montecarlo, needed for stochastic intervention

      # data_cf <- data.table(matrix(0, nrow = 10000, ncol = ncol(data)))
      # setnames(data_cf, colnames(data))
      # data_cf$id = 1:10000
      # data_cf$t = 0
      # task_cf = ltmle3_Task$new(data_cf, tmle_task$npsem, force_at_risk = T)
      # lik = self$cf_likelihood
      # sampler = Sampler$new(lik, c("W",  "A"), "W")
      # tsk = sampler$sample(task_cf, "W", "A")
      # #Make sure task is up to date with targetting
      # tlik$update_task(tsk)
      # liks = lik$get_likelihood(tsk, "processN")
      # liks_surv = liks[, cbind(t,cumprod(.SD)), by = id, .SDcols = "processN"]
      # liks_surv = liks_surv[, lapply(.SD, mean), by = t, .SDcols = c("processN")]
      liks_surv <- 0
      IC =  private$.D_cache[[tmle_task$uuid]]
      if(is.null(IC)){
        self$clever_covariates(tmle_task = tmle_task, fold_number = "full")
      }
      result <- list(psi = liks_surv, IC =  private$.D_cache[[tmle_task$uuid]])
      return(result)
      },
    clever_covariates = function(tmle_task, fold_number = "full"){
      self$clever_covariates_internal(tmle_task, fold_number, subset_times = TRUE)
    },
    get_EIC_component = function(task, node) {
      return(private$.D_cache[[task$uuid]][[node]])
    }

  ),

  active = list(
    # TODO: modify
    name = function() {
      param_form <- sprintf("E[P(T > %s|%s, W)]", self$times, self$cf_likelihood$name)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    times = function() {
      return(private$.times)
    },
    target_times = function() {
      return(private$.target_times)
    }
  ),
  private = list(
    .D_cache = list(),
    .type = "survival",
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .times = NULL,
    .target_times = NULL
  )
)
