
#' Get Inference for Differentiable Functions of Parameters
#'
#' Using the functional delta method, get estimates and inference for a smooth function of parameters
#'
#' @param delta_param \code{list} of functions, defines the parameter
#' @param estimates \code{list}, TMLE estimates of parameter and ICs from \code{\link{tmle3_Fit}$estimates}
#' @param param_names \code{character} vector, the names of the parameters
#' @return \code{list}, delta method estimates of parameter and ICS
#' @export
delta_method <- function(delta_param, estimates, param_names=NULL) {
  psis <- lapply(estimates, `[[`, "psi")
  ICs <- lapply(estimates, `[[`, "IC")
  psi <- delta_param$f(psis)
  IC <- delta_param$df(psis, ICs)
  name <- delta_param$name(param_names)
  list(psi = psi, IC = IC, name = name, transform = delta_param$transform)
}

# Linear Contrast EY1-EY0
f_contrast <- function(x) {
  x[[2]] - x[[1]]
}

df_contrast <- function(x, dx) {
  dx[[2]] - dx[[1]]
}

#' PAR = Linear Contrast EY-EY0
#' @export
delta_param_contrast <- list(
  name = function(names) sprintf("%s - %s", names),
  f = f_contrast,
  df = df_contrast
)

#' PAR = Linear Contrast EY-EY0
#' @export
delta_param_PAR <- list(
  name = function(names) sprintf("PAR(%s)", names[[1]]),
  f = f_contrast,
  df = df_contrast
)
# Risk Ratio EY1/EY0
f_log_rr <- function(x) {
  x[[2]] / x[[1]]
}

df_log_rr <- function(x, dx) {
  dx[[2]] / x[[2]] - dx[[1]] / x[[1]]
}

rr_transform <- exp

#' Risk Ratio EY1/EY0
#' @export
delta_param_RR <- list(
  name = function(names) sprintf("RR(%s/%s)", names[[2]], names[[1]]),
  f = f_log_rr,
  df = df_log_rr,
  transform = rr_transform
)

paf_transform <- function(x) {
  1 - exp(-x)
}

#' PAF = 1 - (1/RR(EY/E0))
#' @export
delta_param_PAF <- list(
  name = function(names) sprintf("PAF(%s)", names[[1]]),
  f = f_log_rr,
  df = df_log_rr,
  transform = paf_transform
)
