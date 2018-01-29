#' @export
delta_method <- function(estimates, f, df){
  psis <- lapply(estimates, `[[`, "psi")
  ICs <- lapply(estimates, `[[`, "IC")
  psi <- f(psis)
  IC <- df(ICs)
  list(psi=psi, IC=IC)
}

#' @export
f_contrast <- function(x){
  x[[2]] - x[[1]]
}

# todo: integrate with tmle_fit methods
#' @export
summary_from_estimate <- function(estimate){
  
  ED2 <- mean(estimate$IC^2)
  se <- sqrt(ED2)/sqrt(length(estimate$IC))
  ci <- wald_ci(estimate$psi, se)
  summary_dt <- data.table(estimate$psi, se, ci)
  setnames(summary_dt, c("tmle_est", "se", "lower","upper"))
  
  return(summary_dt)
}