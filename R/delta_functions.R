# Linear Contrast EY1-EY0
f_contrast <- function(x, dx) {
  x[[2]] - x[[1]]
}

df_contrast <- function(x, dx) {
  dx[[2]] - dx[[1]]
}

#' PAR = Linear Contrast EY1-EY0
#' @export
delta_param_ATE <- list(
  type = "ATE",
  name = function(names) sprintf("%s - %s", names[[2]], names[[1]]),
  f = f_contrast,
  df = df_contrast
)

#' PAR = Linear Contrast EY-EY0
#' @export
delta_param_PAR <- list(
  type = "PAR",
  name = function(names) sprintf("PAR(%s)", names[[1]]),
  f = f_contrast,
  df = df_contrast
)

# Risk Ratio EY1/EY0
f_log_rr <- function(x, dx) {
  log(x[[2]]) - log(x[[1]])
}

df_log_rr <- function(x, dx) {
  dx[[2]] / x[[2]] - dx[[1]] / x[[1]]
}

rr_transform <- exp

#' Risk Ratio EY1/EY0
#' @export
delta_param_RR <- list(
  type = "RR",
  name = function(names) sprintf("RR(%s/%s)", names[[2]], names[[1]]),
  f = f_log_rr,
  df = df_log_rr,
  transform = rr_transform
)

# Odds Ratio odds(Y1)/odds(Y0)
f_log_or <- function(x, dx) {
  log(x[[2]] / (1 - x[[2]])) - log(x[[1]] / (1 - x[[1]]))
}

df_log_or <- function(x, dx) {
  dx[[2]] / (x[[2]] * (1 - x[[2]])) - dx[[1]] / (x[[1]] * (1 - x[[1]]))
}

or_transform <- exp

#' Odds Ratio odds(Y1)/odds(Y0)
#' @export
delta_param_OR <- list(
  type = "OR",
  name = function(names) sprintf("OR(%s/%s)", names[[2]], names[[1]]),
  f = f_log_or,
  df = df_log_or,
  transform = or_transform
)

paf_transform <- function(x) {
  1 - exp(-x)
}

#' PAF = 1 - (1/RR(EY/E0))
#' @export
delta_param_PAF <- list(
  type = "PAF",
  name = function(names) sprintf("PAF(%s)", names[[1]]),
  f = f_log_rr,
  df = df_log_rr,
  transform = paf_transform
)
