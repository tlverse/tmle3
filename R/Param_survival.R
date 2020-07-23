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
      
      private$.cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
      private$.target_times <- target_times
      times <- sort(unique(observed_likelihood$training_task$time))
      private$.times <- times
      if(is.null(target_times)){
        private$.targeted <- rep(TRUE, length(times))
      } else {
        private$.targeted <- times %in% target_times
      }
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
      intervention_nodes <- names(self$intervention_list)
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # I(A=1)
      cf_pA <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      pN <- self$observed_likelihood$get_likelihoods(tmle_task, "N", fold_number)
      # TODO: make bound configurable
      pN <- bound(pN, 0.005)
      
      
      pA_c <- self$observed_likelihood$get_likelihoods(tmle_task, "A_c", fold_number)

      
      
      time <- tmle_task$time
      id <- tmle_task$id
      long_order <- order(id,time)
      
      pA_mat <- self$long_to_mat(pA,id,time)
      t_mat <- self$long_to_mat(time,id,time)
      
      cf_pA_mat <- self$long_to_mat(cf_pA,id,time)
      pN_mat <- self$long_to_mat(pN,id,time)
      pA_c_mat <- self$long_to_mat(pA_c,id,time)
      SN_mat <- self$hm_to_sm(pN_mat)
      SA_c_mat <- self$hm_to_sm(pA_c_mat)
      
      # fix t-1
      SA_c_mat <- cbind(1,SA_c_mat[,-ncol(SA_c_mat)])

      ks <- sort(unique(time))
      
      hk_all <- lapply(ks,function(k){
        Ikt <- k <= t_mat
        SN_mat_k <- matrix(SN_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        SA_c_mat_k <- matrix(SA_c_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        hk <- -1 * ((cf_pA_mat*Ikt)/(pA_mat*SA_c_mat_k))*(SN_mat/SN_mat_k)
      })
      
      # TODO: this might need to be reordered
      HA <- do.call(rbind, hk_all)
      
      if(subset_times & !is.null(self$target_times)){
        HA[,!(ks%in%self$target_times)] = 0
      }
      return(list(N = HA))
    },
    clever_covariates = function(tmle_task, fold_number = "full"){
      self$clever_covariates_internal(tmle_task, fold_number, subset_times = TRUE)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      # TODO: return format
      # TODO: share work between this and the IC code
      HA <- self$clever_covariates_internal(tmle_task, fold_number, subset_times = FALSE)[["N"]]

      time <- tmle_task$time
      id <- tmle_task$id
      
      pN1 <-self$observed_likelihood$get_likelihoods(cf_task, "N", fold_number)
      
      # TODO: make bound configurable
      pN1 <- bound(pN1, 0.005)
      pN1_mat <- self$long_to_mat(pN1,id,time)
      SN1_mat <- self$hm_to_sm(pN1_mat)
      psi <- colMeans(SN1_mat)
      T_tilde <- tmle_task$get_tmle_node("T_tilde")
      Delta <- tmle_task$get_tmle_node("Delta")
      k <- time
      Ittkd <- (T_tilde == k) & (Delta==1)
      Ittk <- (T_tilde >= k)
      
      resid <- as.vector(Ittkd - (Ittk * pN1))
      D1_tk <- HA*resid
      
      # zero out entries that don't contribute to sum
      ts <- sort(unique(k))
      t_mat <- matrix(ts,nrow=nrow(D1_tk),ncol=ncol(D1_tk),byrow = TRUE)
      Itk <- (k<=t_mat)
      D1_tk <- D1_tk * Itk
      D1_tk_dt <- data.table(id=id, k=time, D1_tk)
      
      # sum to IC for 1:N ids
      D1 <- D1_tk_dt[,lapply(.SD,sum),by=list(id),.SDcols=as.character(ts)]
      D1 <- as.matrix(D1[,-1,with=FALSE])
    
      psi_mat <- matrix(psi,nrow=nrow(D1),ncol=ncol(D1),byrow=TRUE)
      D2 <- SN1_mat - psi_mat
      
      IC <- D1+D2
      
      # copy IC to make it match the observation structure
      # TODO: consider if this is the best approach
      IC_id <- sort(unique(id))
      IC_long <- IC[match(id,IC_id),]
      result <- list(psi = psi, IC = IC_long)
      return(result)
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
    .type = "survival",
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .times = NULL,
    .target_times = NULL
  )
)
