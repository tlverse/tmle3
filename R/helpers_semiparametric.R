
get_beta <- function(W, A, formula, Q1, Q0, family, weights = NULL) {
  W <- as.matrix(W)
  if (is.null(weights)) {
    weights <- rep(1, nrow(W))
  }
  V <- model.matrix(formula, as.data.frame(W))
  Q <- ifelse(A == 1, Q1, Q0)
  beta <- suppressWarnings(coef(glm.fit(A * V, Q, offset = family$linkfun(Q0), family = family, intercept = F, weights = weights)))
  return(beta)
}

project_onto_model <- function(W, A, formula, Q1, Q0, family, weights = NULL) {
  beta <- get_beta(W, A, formula, Q1, Q0, family, weights)
  V <- model.matrix(formula, as.data.frame(W))
  Q1 <- family$linkinv(family$linkfun(Q0) + V %*% beta)
  Q <- ifelse(A == 1, Q1, Q0)
  return(cbind(Q0, Q1, Q))
}
