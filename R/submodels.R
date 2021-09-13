
# To port to sl3 at some point:

#' Log likelihood loss for outcomes between 0 and 1
#'
#' @param estimate prediction
#' @param observed observed outcome
#' @export
loss_loglik_binomial <- function(estimate, observed) {
  # loss <- -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
  loss <- -1 * (observed * log(estimate) + (1 - observed) * log(1 - estimate))
  return(loss)
}
#' log likelihood loss
#' @param estimate prediction
#' @param observed observed outcome
#' @export
loss_loglik <- function(estimate, observed) {
  loss <- -1 * log(estimate)
  return(loss)
}

#' Poisson/log-linear loss for nonnegative variables
#'
#' @param estimate prediction
#' @param observed observed outcome
#' @export
loss_poisson <- function(estimate, observed) {
  loss <- estimate - observed * log(estimate)
  return(loss)
}








#' Generate Fluctuation Submodel from \code{family} object.
#'
#' @param family ...
#'
#' @export
#


generate_submodel_from_family <- function(family) {
  linkfun <- family$linkfun
  linkinv <- family$linkinv
  submodel <- function(eps, offset, X, observed) {
    linkinv(linkfun(offset) + X %*% eps)
  }
  return(submodel)
}


#' Logistic Submodel Fluctuation for likelihood (not conditional means)
#'
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#' @importFrom stats plogis qlogis
#'
#' @export
#
submodel_logistic_switch <- function(eps, offset, X, observed) {
  offset <- ifelse(observed == 1, offset, 1 - offset)
  output <- stats::plogis(stats::qlogis(offset) + X %*% eps)
  output <- ifelse(observed == 1, output, 1 - output)
}



#' Generate loss function loss from family object or string
#' @param family ...
#' @export
generate_loss_function_from_family <- function(family) {
  if (!is.character(family)) {
    family <- family$family
  }
  if (!(family %in% c("poisson", "gaussian", "binomial"))) {
    stop("Unsupported family object.")
  }
  if (family == "poisson") {
    return(loss_poisson)
  } else if (family == "gaussian") {
    return(loss_squared_error)
  } else if (family == "binomial") {
    return(loss_loglik_binomial)
  }
}


#' Main maker of submodel specs.
#' @param name ...
#' @export
make_submodel_spec <- function(name, family = NULL, submodel_function = NULL, loss_function = NULL) {
  if (is.null(submodel_function) && inherits(submodel_function, "family")) {
    submodel_function <- generate_submodel_from_family(submodel_function)
  } else if (is.null(submodel_function) && !is.null(family)) {
    submodel_function <- generate_submodel_from_family(family)
  }
  if (is.null(loss_function) && inherits(loss_function, "family")) {
    generate_loss_function_from_family(loss_function)
  } else if (is.null(loss_function) && !is.null(family)) {
    loss_function <- generate_loss_function_from_family(family)
  }
  return(list(name = name, family = family, submodel_function = submodel_function, loss_function = loss_function))
}

#' Main getter for submodel specs.
#' @param name Either a name for a submodel spec obtainable from environment (name -->  get(paste0("submodel_spec_",name))}), a family object or string, or a string of the form "family_link" (e.g. "binomial_logit").
#' @export
get_submodel_spec <- function(name) {
  # If list, assume it is already a spec

  if (is.list(name)) {
    return(name)
  }
  output <- NULL
  tryCatch(
    {
      if (inherits(name, "family")) {
        family <- name
      } else {
        split_names <- unlist(strsplit(name, "_"))
        if (length(split_names) == 2) {
          family <- get(split_names[1])(link = split_names[2])
        } else {
          family <- get(split_names[1])()
        }
      }
      output <- make_submodel_spec(name, family)
    },
    error = function(...) {
      print(...)
      try({
        output <<- get(paste0("submodel_spec_", name))
      })
    }
  )
  if (is.null(output)) {
    stop(paste0("Argument name was not a valid family nor was `submodel_spec_", name, "` found in the environment."))
  }
  return(output)
}

#' Submodel for binary outcomes where "initial" is a likelihood and not a conditional mean (e.g. for Param_ATC and Param_ATT for updating node `A`).
#' @export
submodel_spec_logistic_switch <- list(name = "logistic_switch", family = function() {
  stop("Does not support family-based updating. Please use optim instead.")
}, submodel_function = submodel_logistic_switch, loss_function = loss_loglik)
