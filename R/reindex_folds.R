reindex <- function(index, subset) {
  matches <- match(index, subset)
  reindexed <- matches[!is.na(matches)]
  return(reindexed)
}

subset_fold <- function(fold, subset) {
  origami::make_fold(
    origami::fold_index(),
    reindex(origami::training(), subset),
    reindex(origami::validation(), subset)
  )
}

subset_folds <- function(folds, subset) {
  subsetted <- lapply(folds, subset_fold, subset)
  return(subsetted)
}
