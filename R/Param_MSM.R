#' Stratified Parameter Estimates via MSM
#'
#' @section Current Issues:
#' \itemize{
#'   \item clever covariates doesn't support updates; always uses initial (necessary for iterative TMLE, e.g. stochastic intervention)
#'   \item clever covariate gets recalculated all the time (inefficient)
#' }
#' @importFrom R6 R6Class
#' @importFrom stats glm predict
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importFrom stringr str_split str_replace_all
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_MSM, observed_likelihood, strata_variable, ...)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{msm}}{form of the MSM. Default is "A + V", consistent with the default of 
#'     \code{intervention_node} and \code{strata_name}.
#'     }
#'     \item{\code{weight}}{"Cond.Prob.", "Unif." or custom input function. 
#'     Note that custom function should support vector input. Default is "Cond.Prob.".
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{covariate_node}}{character, the name of the node that should be treated as the covariate
#'     }
#'     \item{\code{intervention_node}}{character, the name of the node that should be treated as the intervention
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this intervention
#'     }
#' }
#' @export
Param_MSM <- R6Class(
  classname = "Param_MSM",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, strata_variable, strata_name = "V",
                          msm = "A + V", weight = "Cond.Prob.", weight_ub = 1/0.025,
                          continuous_intervention = FALSE, intervention_values = NULL, n_samples = 30,
                          ...,
                          covariate_node = "W", intervention_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.strata_variable <- strata_variable
      private$.strata_name <- strata_name
      private$.weight_ub <- weight_ub
      private$.continuous_intervention <- continuous_intervention
      private$.covariate_node <- covariate_node
      private$.intervention_node <- intervention_node
      
      if (continuous_intervention) {
        private$.n_samples <- n_samples
        #private$.U <- runif(n_samples)
        private$.U <- seq(0, 1, length.out = self$n_samples)
        private$.intervention_values <- setNames(self$intervention_node, self$intervention_node) # not needed but stored for simplicity
      } else {
        if (is.null(intervention_values)) {
          intervention_values <- observed_likelihood$factor_list[[self$intervention_node]]$variable_type$levels
        }
        private$.intervention_values <- setNames(intervention_values, paste0(self$intervention_node, "_", intervention_values))
        # cf_likelihoods (not needed for internal tasks)
        private$.cf_likelihoods <- lapply(self$intervention_values, function(a) {
          intervention_list <- define_lf(LF_static, self$intervention_node, value = a)
          cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
          return(cf_likelihood)
        })
        names(private$.cf_likelihoods) <- names(self$intervention_values)
      }
      
      # terms of MSM
      msm <- tail(str_split(msm, "[ ]*~[ ]*")[[1]], n = 1)
      private$.msm_terms <- str_split(msm, "[ ]*\\+[ ]*")[[1]]
      
      private$.msm_terms_all <- c()
      for (t in self$msm_terms) {
        if (grepl(self$intervention_node, t)) {
          private$.msm_terms_all <- c(private$.msm_terms_all, str_replace_all(t, self$intervention_node, names(self$intervention_values)))
        } else {
          private$.msm_terms_all <- c(private$.msm_terms_all, t)
        }
      }
      
      # process h(A, V)
      if (weight == "Cond.Prob.") {
        # current: P(A|V,W)
        # TODO: refit P(A|V)
        private$.weight <- function(tmle_task, fold_number = "full") {
          h = self$observed_likelihood$get_likelihoods(tmle_task, self$intervention_node, fold_number)
          if (!is.null(self$weight_ub)) {
            h = pmin(h, self$weight_ub)
          }
          h
        }
      } else if (weight == "Unif.") {
        private$.weight <- function(tmle_task, fold_number = "full") {
          h = rep(1, tmle_task$nrow)
          if (!is.null(self$weight_ub)) {
            h = pmin(h, self$weight_ub)
          }
          h
        }
      } else {
        if (!is.function(weight)) {
          stop("weight should be either a valid string or a function h(A, V). \n")
        }
        private$.weight <- function(tmle_task, fold_number = "full") {
          # TODO: fold?
          A <- tmle_task$get_tmle_node(self$intervention_node)
          V <- tmle_task$get_tmle_node(self$covariate_node)[[self$strata_variable]]
          h = sapply(1:length(A), function(i) weight(A[i], V[i]))
          if (!is.null(self$weight_ub)) {
            h = pmin(h, self$weight_ub)
          }
          h
        }
      }
    },
    
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      g <- self$observed_likelihood$get_likelihoods(tmle_task, self$intervention_node, fold_number)
      h <- private$.weight(tmle_task, fold_number)
      if (self$continuous_intervention) {
        A <- as.matrix(tmle_task$get_tmle_node(self$intervention_node))
      } else {
        A <- self$get_intervention_indicators(tmle_task)
      }
      V <- as.matrix(tmle_task$get_data(, self$strata_variable))
      
      msm_terms <- self$msm_terms %>% str_replace_all(self$intervention_node, "A") %>% str_replace_all(self$strata_name, "V")
      phi <- do.call(cbind, lapply(msm_terms, function(t) eval(parse(text = t))))
      
      H1 <- c(h / g) * phi
      
      colnames(H1) <- self$msm_terms_all
      return(list(Y = H1))
    },
    
    estimates = function(tmle_task = NULL, fold_number = "full") {
      # TODO: estimate variance of MC integration
      n_obs <- tmle_task$nrow
      
      Y <- as.matrix(tmle_task$get_tmle_node(self$outcome_node))
      Q <- self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node, fold_number)
      V <- as.vector(as.matrix(tmle_task$get_data(, self$strata_variable)))
      H1 <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
      
      if (self$continuous_intervention) {
        A_range <- self$observed_likelihood$factor_list[[self$intervention_node]]$learner$get_outcome_range(tmle_task$get_regression_task(self$outcome_node), fold_number)
        A_vals <- t(apply(A_range, 1, function(r) r[1] + (r[2] - r[1]) * self$U))
        # Generate counterfactual tasks for each sample of A:
        cf_tasks <- apply(A_vals, 2, function(a) {
          newdata <- data.table(A = a)
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
          return(cf_task)
        })
      } else {
        A_vals <- diag(length(self$intervention_values))
        # Generate counterfactual tasks for each value of A:
        cf_tasks <- lapply(self$cf_likelihoods, function(cf_likelihood) {
          cf_task <- cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
          return(cf_task)
        })
      }
      
      n_treats <- ncol(A_vals)
      
      QA <- sapply(cf_tasks, self$observed_likelihood$get_likelihood, self$outcome_node, fold_number)
      hA <- sapply(cf_tasks, private$.weight, fold_number)
      
      # psi
      Q_ext <- matrix(QA, ncol = 1, byrow = FALSE)
      if (self$continuous_intervention) {
        V_ext <- rep(V, times = self$n_samples)
        A_ext <- matrix(A_vals, ncol = 1, byrow = FALSE)
      } else {
        V_ext <- rep(V, times = length(self$intervention_values))
        A_ext <- sapply(data.frame(A_vals), rep, each = n_obs)
      }
      h_ext <- matrix(hA, ncol = 1, byrow = FALSE)
      
      msm_terms <- self$msm_terms %>% str_replace_all(self$intervention_node, "A_ext") %>% str_replace_all(self$strata_name, "V_ext")
      phi_ext <- do.call(cbind, lapply(msm_terms, function(t) eval(parse(text = t))))
      
      regress_table <- data.table(cbind(Q_ext, phi_ext))
      colnames(regress_table) <- c("Q", colnames(H1))
      formula <- as.formula(paste0("Q~", 
                                   paste(colnames(H1), collapse = "+"),
                                   "-1"))
      if (length(levels(factor(Y))) == 2 && all.equal(levels(factor(Y)), c("0", "1"))) {
        family = "binomial"
      } else {
        family = "gaussian"
      }
      suppressWarnings(model <- glm(formula, family = family, data = regress_table, weights = h_ext))
      
      psi <- model$coefficients
      names(psi) <- colnames(H1)
      
      # ic
      weighted_res <- c(residuals(model, type = "response") * h_ext)
      
      H2A <- split(as.data.table(weighted_res * phi_ext), 
                       factor(rep(1:n_treats, each = n_obs)))
      H2 <- Reduce('+', lapply(H2A, as.matrix))
      
      IC <- H1 * c(Y - Q) + H2
      
      # normalization
      if (!any(is.na(IC))) {
        if (family == "binomial") {
          m_ext <- matrix(predict(model, type = "response"))
          dm_ext <- m_ext * (1 - m_ext)
        } else {
          dm_ext <- matrix(1, nrow = n_obs * n_treats, ncol = 1)
        }
        
        M <- t(c(dm_ext * h_ext) * phi_ext) %*% phi_ext / n_obs
        
        Minv <- try(solve(M))
        if (identical(class(Minv), "try-error")) {
          warning("Inference unavailable: normalizing matrix not invertible. IC not normalized. \n")
        } else {
          IC <- IC %*% Minv
        }
      }
      
      colnames(IC) <- colnames(H1)
      result <- list(psi = psi, IC = IC)
      return(result)
    },
    
    get_intervention_indicators = function(tmle_task, intervention_node = self$intervention_node) {
      intervention_variable <- tmle_task$npsem[[intervention_node]]$variables
      A <- tmle_task$get_data(, intervention_variable)
      intervention <- factor(A[[intervention_variable]])
      
      combined <- A[, "indicators":= 1][, "level":= intervention]
      intervention_indicators <- t(sapply(combined$level, function(a) as.numeric(a == self$intervention_values)))
      colnames(intervention_indicators) <- names(self$intervention_values)
      return(intervention_indicators)
    }
  ),
  
  active = list(
    name = function() {
      param_form <- sprintf("beta[%s]", self$msm_terms_all)
      return(param_form)
    },
    continuous_intervention = function() {
      return(private$.continuous_intervention)
    },
    weight_ub = function() {
      return(private$.weight_ub)
    },
    intervention_values = function() {
      return(private$.intervention_values)
    },
    n_samples = function() {
      return(private$.n_samples)
    },
    strata_variable = function() {
      return(private$.strata_variable)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    strata_name = function() {
      return(private$.strata_name)
    },
    covariate_node = function() {
      return(private$.covariate_node)
    },
    intervention_node = function() {
      return(private$.intervention_node)
    },
    U = function() {
      return(private$.U)
    },
    cf_likelihoods = function() {
      return(private$.cf_likelihoods)
    },
    msm_terms = function() {
      return(private$.msm_terms)
    },
    msm_terms_all = function() {
      return(private$.msm_terms_all)
    }
  ),
  
  private = list(
    .type = "stratified MSM",
    .continuous_intervention = NULL,
    .weight_ub = NULL,
    .intervention_values = NULL,
    .n_samples = NULL,
    .strata_variable = NULL,
    .strata_name = NULL,
    .weight = NULL,
    .covariate_node = NULL,
    .intervention_node = NULL,
    .U = NULL,
    .cf_likelihoods = NULL,
    .msm_terms = NULL,
    .msm_terms_all = NULL
  )
)