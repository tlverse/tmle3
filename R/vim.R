utils::globalVariables(c(
  "Z_stat", "tmle_est", "se", "p_nz", "p_nz_corrected",
  "A"
))

#' @importFrom stats p.adjust pnorm
#' @importFrom foreach foreach "%do%"
#' @export
tmle3_vim <- function(tmle_spec, data, node_list, learner_list = NULL, adjust_for_other_A = TRUE) {
  A_nodes <- node_list$A

  vim_A <- A_nodes[[1]]
  all_estimates <- foreach::foreach(vim_A = A_nodes) %do% {
    vim_node_list <- node_list
    if (adjust_for_other_A) {
      vim_node_list$W <- unique(c(vim_node_list$W, vim_node_list$A))
    }

    vim_node_list$W <- setdiff(vim_node_list$W, vim_A)
    vim_node_list$A <- vim_A

    # TODO: consider if there should be a shared outcome model here
    vim_fit <- tmle3(tmle_spec, data, vim_node_list, learner_list)


    # pull off final estimate
    # TODO: allow estimate to be user-selectable
    estimate <- vim_fit$summary[length(vim_fit$estimates)]
    estimate$A <- vim_A
    estimate$W <- list(list(vim_node_list$W))
    estimate
  }

  vim_results <- rbindlist(all_estimates)
  vim_results[, Z_stat := tmle_est / se]

  # TODO: think about VIM parameters where H0 is something other than 0
  vim_results[, p_nz := stats::pnorm(-1 * abs(Z_stat))]
  vim_results[, p_nz_corrected := p.adjust(p_nz, method = "BH")]
  vim_results <- vim_results[order(p_nz_corrected)]
  vim_results <- vim_results[, A := factor(A, levels = A)]
  return(vim_results)
}

#' @import ggplot2
#' @export
plot_vim <- function(vim_results) {
  ggplot(
    vim_results,
    aes_string(
      x = "psi_transformed", xmin = "lower_transformed",
      xmax = "upper_transformed", y = "A"
    )
  ) +
    geom_point() +
    geom_errorbarh() +
    theme_bw() +
    xlab("Importance Measure") + ylab("Variable") +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_y_discrete(limits = rev(levels(vim_results$A)))
}
