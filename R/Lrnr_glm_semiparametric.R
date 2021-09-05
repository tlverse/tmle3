#' Semiparametric Generalized Linear Models
#'
#' This learner provides fitting procedures for generalized linear models using
#' \code{\link[stats]{glm.fit}}.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom stats glm predict family
#'
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{...}}{Parameters passed to \code{\link[stats]{glm}}.}
#' }
#'
#
Lrnr_glm_semiparametric <- R6Class(
  classname = "Lrnr_glm_semiparametric", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(formula_sp, lrnr_baseline, interaction_variable = "A", family = NULL, append_interaction_matrix = TRUE, return_matrix_predictions = F, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "weights", "offset"),

    .train = function(task) {

      args <- self$params
      append_interaction_matrix <- args$append_interaction_matrix
      outcome_type <- self$get_outcome_type(task)
      trt <- args$interaction_variable
      if(is.null(trt)) {
        A <- rep(1, task$nrow)
      } else {
        A <- unlist(task$get_data(,trt))
      }
      if(!all(A %in% c(0,1)) && !is.null(trt)) {
        binary <- FALSE
      } else {
        binary <- TRUE
      }
      family <- args$family
      lrnr_baseline <- args$lrnr_baseline
      formula <- args$formula_sp
      if (is.null(family)) {
        family <- outcome_type$glm_family(return_object = TRUE)
      }
      # Interaction design matrix
      Y <- task$Y
      V <- model.matrix(formula, task$data)
      colnames(V) <- paste0("V", 1:ncol(V))

      covariates <- setdiff(task$nodes$covariates, trt)

      if(!append_interaction_matrix && binary) {
        task_baseline <- task$next_in_chain(covariates = covariates)
        lrnr_baseline <- lrnr_baseline$train(task_baseline[A==0])
        Q0 <- lrnr_baseline$predict(task_baseline)
        beta <- coef(glm.fit(A*V, Y, offset = family$linkfun(Q0), intercept = F, weights = task$weights, family = family))
        Q1 <- family$linkinv(family$linkfun(Q0) + V%*%beta)
        Q <- ifelse(A==1, Q1, Q0)
      } else {

        covariates <- setdiff(task$nodes$covariates, trt)

        if(append_interaction_matrix) {
          AV <- as.data.table(A*V)
          X <- cbind(task$X[,covariates, with = F], AV)
          X0 <- cbind(task$X[,covariates, with = F], 0*V)
        } else {
          X <- cbind(task$X[,covariates, with = F], A)
          X0 <- cbind(task$X[,covariates, with = F], A*0)
        }


        column_names <- task$add_columns(X)
        task_baseline <- task$next_in_chain(covariates = colnames(X), column_names = column_names )

        column_names <- task$add_columns(X0)
        task_baseline0 <- task$next_in_chain(covariates = colnames(X0), column_names = column_names )

        lrnr_baseline <- lrnr_baseline$train(task_baseline)
        Q <- lrnr_baseline$predict(task_baseline)
        Q0 <- lrnr_baseline$predict(task_baseline0)
        # Project onto model

        beta <- coef(glm.fit(A*V, Q, offset = family$linkfun(Q0), intercept = F, weights = task$weights, family = family))

      }

       fit_object = list(beta = beta, lrnr_baseline = lrnr_baseline, covariates = covariates, family = family, formula = formula,
                         append_interaction_matrix = append_interaction_matrix, binary = binary, task_baseline = task_baseline)
      return(fit_object)
    },
    .predict = function(task) {
      fit_object <- self$fit_object
      append_interaction_matrix <- fit_object$append_interaction_matrix
      binary <- fit_object$binary
      beta <- fit_object$beta
      lrnr_baseline <- fit_object$lrnr_baseline
      covariates <- fit_object$covariates
      family <- fit_object$family
      formula <- fit_object$formula

      trt <- self$params$interaction_variable
      if(is.null(trt)) {
        A <- rep(1, task$nrow)
      } else {
        A <- unlist(task$get_data(,trt))
      }
      V <- model.matrix(formula, task$data)
      colnames(V) <- paste0("V", 1:ncol(V))


      if(!append_interaction_matrix && binary) {
        task_baseline <- task$next_in_chain(covariates = covariates)
        Q0 <- lrnr_baseline$predict(task_baseline)
      } else {
        if(append_interaction_matrix) {
          X0 <- cbind(task$X[,covariates, with = F], 0*V)
        } else {
          X0 <- cbind(task$X[,covariates, with = F], 0)
        }
        column_names <- task$add_columns(X0)
        task_baseline0 <- task$next_in_chain(covariates = colnames(X0), column_names = column_names )
        Q0 <-  lrnr_baseline$predict(task_baseline0)
      }

      Q1 <- family$linkinv(family$linkfun(Q0) + V%*%beta)
      Q <- family$linkinv(family$linkfun(Q0) + A*V%*%beta)
      if(self$params$return_matrix_predictions && binary) {
        predictions <- cbind(Q0,Q1,Q)
        colnames(predictions) <- c("A=0", "A=1", "A")
        predictions <- sl3::pack_predictions(cbind(Q0,Q1))
      } else {
        predictions <- Q
      }


      return(predictions)
    }
  )
)
