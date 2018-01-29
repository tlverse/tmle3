wald_ci <- function(est, se, alpha=0.95){
  z <- qnorm((1+alpha)/2)
  lower <- est - z * se
  upper <- est + z * se
  ci <- cbind(lower, upper)
  return(ci)  
}

#' Fit Targeted Maximum Likelihood Estimator
#'
#' @param likelihood ...
#' @param task ...
#' @param param ...
#' @param lrnr_submodel ...
#'
#' @export
#' @import sl3
#' @import data.table
tmle3_Fit <- R6Class(
  classname = "tmle3_Fit",
  public = list(
    initialize = function(tmle_task, likelihood, tmle_params, updater, max_it=100){
      private$.tmle_task <- tmle_task
      private$.likelihood <- likelihood
      private$.tmle_params <- tmle_params
      private$.updater <- updater
      initial_psi <- sapply(self$tmle_params, function(tmle_param) tmle_param$estimates(self$tmle_task)$psi)
      private$.initial_psi <- initial_psi
      private$.tmle_fit(max_it)
    },
    print = function(){
      cat(sprintf("A tmle3_Fit that took %s step(s)\n",self$steps))
      print(self$summary)
    }
  ),
  active = list(
    tmle_task = function() {
      return(private$.tmle_task)
    },
    likelihood = function() {
      return(private$.likelihood)
    },
    tmle_params = function() {
      return(private$.tmle_params)
    },
    tmle_param_names = function() {
      if(is.null(private$.tmle_param_names)){
        private$.tmle_param_names <- sapply(self$tmle_params, `[[`, "name")  
      }
      return(private$.tmle_param_names)
    },
    updater = function() {
      return(private$.updater)
    },
    steps = function(){
      return(private$.steps)
    },
    ED = function(){
      ED <- private$.ED
      # names(ED) <- self$tmle_param_names
      return(ED)
    },
    initial_psi = function(){
      initial_psi <- private$.initial_psi
      # names(initial_psi) <- self$tmle_param_names
      return(initial_psi)
    },
    estimates = function(){
      estimates <- private$.estimates
      # names(estimates) <- self$tmle_param_names
      return(estimates)
    },
    psi = function(){
      psi <- sapply(self$estimates, `[[`, "psi")
      # names(psi) <- self$tmle_param_names
      return(psi)
    },
    IC = function(){
      IC <- sapply(self$estimates, `[[`, "IC")
      # names(IC) <- self$tmle_param_names
      return(IC)
    },
    se = function(){
      IC <- self$IC
      ED2 <- apply(IC, 2, function(x)mean(x^2))
      se <- sqrt(ED2)/sqrt(self$tmle_task$nrow)
      return(se)
    },
    ci = function(){
      return(wald_ci(self$psi, self$se))
    },
    summary = function(){
      summary_dt <- as.data.table(list(self$tmle_param_names, self$initial_psi, self$psi, self$se, self$ci))
      setnames(summary_dt, c("param", "init_est", "tmle_est", "se", "lower","upper"))
      return(summary_dt)
    }
  ),
  private = list(
    .tmle_task = NULL,
    .likelihood = NULL,
    .tmle_params = NULL,
    .tmle_param_names = NULL,
    .updater = NULL,
    .steps = NULL,
    .ED = NULL,
    .initial_psi = NULL,
    .estimates = NULL,
    .tmle_fit = function(max_it=100){
      ED_criterion <- 1 / self$tmle_task$nrow
      
      for (steps in 1:max_it) {
        self$updater$update_step(self$tmle_task, self$likelihood)
        
        estimates <- lapply(self$tmle_params, function(tmle_param) tmle_param$estimates(self$tmle_task))
        ICs <- sapply(estimates, `[[`, "IC")
        ED <- colMeans(ICs)
        if (max(abs(ED)) < ED_criterion) {
          break
        }
      }
      
      private$.ED=ED
      private$.steps=steps
      private$.estimates = estimates
      
    }
  )
)

#' @export
fit_tmle3 <- function(tmle_task, likelihood, tmle_params, updater) {
  tmle3_Fit$new(tmle_task, likelihood, tmle_params, updater)
}

#' @export
tmle3 <- function(tmle_spec, data, node_list, learner_list = NULL) {
  tmle_spec$tmle3(data, node_list, learner_list)
}

#' @export
plot.tmle3_Fit = function(x){
  summary <- x$summary
  est_names <- c("init_est", "tmle_est")
  est_labels <- c("Initial", "TMLE")
  long <- melt(summary, id=c("param", "lower","upper"), measure=est_names)
  long$variable <- est_labels[match(long$variable, est_names)]
  long[variable=="Initial", lower:=NA]
  long[variable=="Initial", upper:=NA]
  ggplot(long, aes(y=param, x=value, xmin=lower, xmax=upper, color=variable))+
    geom_point()+geom_errorbarh(data=long[!is.na(lower)])+theme_bw()+xlab("Value")+ylab("Parameter")+scale_color_discrete("")
}