#' @import ggplot2
#' @export
plot.tmle3_Fit <- function(x, ...) {
  # kludge for Rcmd::check with data.table:
  # see https://github.com/Rdatatable/data.table/issues/850
  variable <- lower <- upper <- NULL

  summary <- x$summary
  est_names <- c("init_est", "tmle_est")
  est_labels <- c("Initial", "TMLE")
  long <- melt(summary, id = c("param", "lower", "upper"), measure = est_names)
  long$variable <- est_labels[match(long$variable, est_names)]
  long[variable == "Initial", lower := NA]
  long[variable == "Initial", upper := NA]
  ggplot(long, aes_(y = ~param, x = ~value, xmin = ~lower, xmax = ~upper, color = ~variable)) +
    geom_point() +
    geom_errorbarh(data = long[!is.na(lower)]) +
    theme_bw() +
    xlab("Value") +
    ylab("Parameter") +
    scale_color_discrete("")
}
