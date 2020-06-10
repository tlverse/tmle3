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
    initialize = function(observed_likelihood, intervention_list, ..., outcome_node) {
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
    },
    reshape_long_data = function(long_data, t_max) {
      n <- length(long_data) / t_max
      # TODO: assume long_data is a list
      rs <- list()
      for (i in 1:t_max) {
        current <- long_data[seq(1 + (i - 1) * n, i * n)]
        rs <- c(rs, list(current))
      }
      rs <- do.call(cbind, rs)
      return(rs)
    },
    hazards_to_survival = function(p_hazards, t_max) {
      # TODO: whether change the input data format to consecutive times for each observation
      n <- length(p_hazards) / t_max
      p_surv <- copy(p_hazards)
      for (i in 1:n) {
        # TODO: make hazard at starting time 0
        temp <- p_hazards[i + n * seq(0, t_max - 1)]
        temp <- c(0, temp)
        p_surv[i + n * seq(0, t_max - 1)] <- cumprod(1 - temp[-length(temp)])
      }
      return(p_surv)
    },
    get_pSN_at_time = function(pS_N, time, t_max, long_format = TRUE) {
      n <- length(pS_N) / t_max
      pS_N_time <- pS_N[seq(1 + (time - 1) * n, time * n)]
      if (long_format) {
        pS_N_time <- rep(pS_N_time, t_max)
      } 
      return(pS_N_time)
    },
    get_single_time_ht = function(time, pA, cf_pA, pS_N, pS_A_c, k_list, t_max) {
      I_k <- ifelse(k_list <= time, 1, 0)
      pS_N_time <- self$get_pSN_at_time(pS_N, time, t_max)
      ht <- -((cf_pA * I_k) / (pA * pS_A_c)) * (pS_N_time / pS_N)
      return(ht)
    },
    get_hazard = function(tmle_task, node, fold_number){
      oH <- tmle_task$get_tmle_node(node)
      pH <- self$observed_likelihood$get_likelihoods(tmle_task, node, fold_number)
      return(ifelse(oH==1,pH,1-pH))
    },
    long_to_mat = function(x,id, time){
      dt <- data.table(id=id,time=time,x=as.vector(x))
      wide <- dcast(dt, id~time, value.var="x")
      mat <- as.matrix(wide[,-1,with=FALSE])
      return(mat)
    },
    hm_to_sm = function(hm){
      t(apply(1-hm,1,cumprod))
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- names(self$intervention_list)
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # I(A=1)
      cf_pA <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      # TODO: whether modify LF_fit_hazards get_density without 1 - preds
      
      pN <- self$get_hazard(tmle_task, "N", fold_number)
      pA_c <- self$get_hazard(tmle_task, "A_c", fold_number)

      t_max <- max(tmle_task$get_tmle_node("T_tilde"))
      
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

      ks <- sort(unique(time))
      
      hk_all <- lapply(ks,function(k){
        Ikt <- k <= t_mat
        SN_mat_k <- matrix(SN_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        SA_c_mat_k <- matrix(SA_c_mat[,k],nrow=nrow(t_mat),ncol=ncol(t_mat))
        hk <- -1 * ((cf_pA_mat*Ikt)/(pA_mat*SA_c_mat_k))*(SN_mat/SN_mat_k)
      })
      
      # TODO: this might need to be reordered
      HA <- do.call(rbind, hk_all)
      
      
      return(list(N = HA))
    },
    get_psi = function(pS_N1, t_max) {
      n <- length(pS_N1) / t_max
      psi <- lapply(seq(t_max), function(t) {
        mean(pS_N1[seq(1 + (t - 1) * n, t * n)])
      })
      return(unlist(psi))
    },
    get_single_time_Dt = function(time, HA, pN1, pS_N1, psi, T_tilde_data_short, Delta_data_short, t_max) {
      n <- length(pS_N1) / t_max
      # TODO: initialize properly
      Dt <- NULL
      for (k in 1:time) {
        # TODO: optimize I1, I2 creation
        I1 <- ifelse(T_tilde_data_short == k & Delta_data_short == 1, 1, 0)
        I2 <- ifelse(T_tilde_data_short >= k, 1, 0)
        # TODO: check
        temp = HA[seq(1 + (k - 1) * n, k * n), time] * 
        (I1 - I2 * pN1[seq(1 + (k - 1) * n, k * n)])
        if (is.null(Dt)) {
          Dt = temp
        } else {
          Dt = Dt + temp
        }
      }
      Dt = Dt + self$get_pSN_at_time(pS_N1, time, t_max, long_format = FALSE) - psi[time]
      return(Dt)
    },
    # get_all_Dt = function(HA, pN1, pS_N1, psi, T_tilde_data_short, Delta_data_short, t_max) {
    #   n <- length(pS_N1) / t_max

    #   all_Dt_mat <- matrix(rep(0, n * t_max), nrow = n, ncol = t_max, byrow = TRUE)
    #   all_Dt <- as.data.frame(all_Dt_mat)
    #   cum_sum1 <- NULL
    #   for (t in 1:t_max) { 
    #     I1 <- ifelse(T_tilde_data_short == t & Delta_data_short == 1, 1, 0)
    #     I2 <- ifelse(T_tilde_data_short >= t, 1, 0)
    #     # TODO: not work because ht differs for different t
    #     part1 <- HA[seq(1 + (t - 1) * n, t * n), t] * 
    #     (I1 - I2 * pN1[seq(1 + (t - 1) * n, t * n)])
    #     part2 <- self$get_pSN_at_time(pS_N1, t, t_max, long_format = FALSE) - psi[t]
    #     if (is.null(cum_sum1)) {
    #       cum_sum1 <- part1
    #     } else {
    #       cum_sum1 = cum_sum1 + part1
    #     }
    #     all_Dt[, t] <- cum_sum1 + part2
    #   }
    #   return(all_Dt)
    # },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      # TODO: return format
      # TODO: share work between this and the IC code
      HA <- self$clever_covariates(tmle_task, fold_number)[["N"]]

      time <- tmle_task$time
      id <- tmle_task$id
      
      pN1 <- pN <- self$get_hazard(cf_task, "N", fold_number)
      pN1_mat <- self$long_to_mat(pN1,id,time)
      SN1_mat <- self$hm_to_sm(pN1_mat)
      psi <- colMeans(SN1_mat)
      T_tilde <- tmle_task$get_tmle_node("T_tilde")
      Delta <- tmle_task$get_tmle_node("Delta")
      k <- time
      Itkd <- (T_tilde == k) & (Delta==1)
      Itk <- (T_tilde >= k)
      resid <- as.vector(Itkd - (Itk * pN1))
      D1_tk <- HA*resid
      D1_tk_dt <- data.table(id=id, k=time, D1_tk)
      # for each id, t, sum D1_tk for k<=t
      ts <- sort(unique(k))
      D1_all <- lapply(ts, function(t){
        col_name <-  names(D1_tk_dt)[t+2]
        temp <- D1_tk_dt[k<=t,sum(.SD),by=list(id),.SDcols=col_name]
        unlist(temp$V1)
      })
      D1 <- do.call(cbind,D1_all)
      psi_mat <- matrix(psi,nrow=nrow(D1),ncol=ncol(D1),byrow=TRUE)
      D2 <- SN1_mat - psi_mat
      
      IC <- D1+D2
      
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    # TODO: modify
    name = function() {
      param_form <- sprintf("E[%s_{%s}]", self$outcome_node, self$cf_likelihood$name)
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
    }
  ),
  private = list(
    .type = "survival",
    .cf_likelihood = NULL
  )
)
