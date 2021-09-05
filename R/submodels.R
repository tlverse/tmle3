#' Generate Fluctuation Submodel from \code{family} object.
#'
#' @param family ...
#'
#' @export
#
generate_submodel_from_family <- function(family) {
  linkfun <- family$linkfun
  linkinv <- family$linkinv
  submodel <- function(eps, offset, X) {
    linkinv(linkfun(offset) + X %*% eps)
  }
  return(submodel)
}


#' Logistic Submodel Fluctuation
#'
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#' @importFrom stats plogis qlogis
#'
#' @export
#
submodel_logistic <- generate_submodel_from_family(binomial())

#' Log likelihood loss for binary variables
#'
#' @param estimate ...
#' @param observed ...
#' @param weights ...
#' @param v ...
#' @export
loss_function_loglik_binomial = function(estimate, observed, weights = NULL, likelihood = NULL, tmle_task = NULL, fold_number = NULL) {
  loss <- -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
  if(!is.null(weights)) {
    loss <- weights * loss
  }
  return(loss)
}

#' Linear (gaussian) Submodel Fluctuation
#'
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#'
#' @export
#
submodel_linear <- generate_submodel_from_family(gaussian())
#' Least-squares loss for binary variables
#'
#' @param estimate ...
#' @param observed ...
#' @param weights ...
#' @param likelihood ...
#' @export
loss_function_least_squares = function(estimate, observed, weights = NULL,  likelihood = NULL, tmle_task = NULL, fold_number = NULL) {
  loss <- (observed - estimate)^2
  if(!is.null(weights)) {
    loss <- weights * loss
  }
  return(loss)
}


#' Log-linear (Poisson) Submodel Fluctuation
#'
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#'
#' @export
#
submodel_exp  <- generate_submodel_from_family(poisson())

#' Poisson/log-linear loss for nonnegative variables
#'
#' @param estimate ...
#' @param observed ...
#' @param weights ...
#' @param likelihood ...
#' @export
loss_function_poisson = function(estimate, observed, weights = NULL, likelihood = NULL, tmle_task = NULL, fold_number = NULL) {
  loss <-  estimate - observed * log(estimate)
  if(!is.null(weights)) {
    loss <- weights * loss
  }
  return(loss)
}

#' Generate loss function loss from family object or string
#' @param family ...
#' @export
generate_loss_function_from_family <- function(family) {
  if(!is.character(family)) {
    family <- family$family
  }
  if(!(family %in% c("poisson", "gaussian", "binomial"))){
    stop("Unsupported family object.")
  }
  if(family == "poisson"){
    return(loss_function_poisson)
  } else if(family == "gaussian"){
    return(loss_function_least_squares)
  } else if(family == "binomial"){
    return(loss_function_loglik_binomial)
  }
}


#' Main maker of submodel specs.
#' @param name ...
#' @export
make_submodel_spec <- function(name, family = NULL, submodel_function  = NULL, loss_function  = NULL) {
  if(is.null(submodel_function)  && inherits(submodel_function, "family")) {
    submodel_function <- generate_submodel_from_family(submodel_function)
  } else if(is.null(submodel_function) && !is.null(family)) {
    submodel_function <- generate_submodel_from_family(family)
  }
  if(is.null(loss_function)  && inherits(loss_function, "family")) {
    generate_loss_function_from_family(loss_function)
  } else if(is.null(loss_function) && !is.null(family)) {
    loss_function <- generate_loss_function_from_family(family)
  }
  return(list(name = name, family = family, submodel_function = submodel_function, loss_function = loss_function))
}

#' Main getter for submodel specs.
#' @param name Either a name for a submodel spec obtainable from environment (name -->  get(paste0("submodel_spec_",name))}), a family object or string, or a string of the form "family_link" (e.g. "binomial_logit").
#' @export
get_submodel_spec <- function(name) {
  output <- NULL
  tryCatch({
    if(inherits(name, "family")) {
      family <- name
    } else {
      split_names <- unlist(strsplit(name, "_"))
      if(length(split_names)==2) {
        family <- get(split_names[1])(link = split_names[2])
      } else {
        family <- get(split_names[1])()
      }
    }
    output <- make_submodel_spec(name, family)
  }, error = function(...) {
    try({output <<- get(paste0("submodel_spec_",name))})
  })
  if(is.null(output)) {
    stop(paste0("Argument name was not a valid family nor was `submodel_spec_", name, "` found in the environment."))
  }
  return(output)
}

