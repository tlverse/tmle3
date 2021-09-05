

submodel_function <- function(eps, offset, X, observed) {
  offset_unpacked <- sl3::unpack_predictions(offset)
  Q0 <- offset_unpacked[[1]]
  Q1 <- offset_unpacked[[2]]
  Q <- offset_unpacked[[3]]

}
