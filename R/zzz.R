#' @importFrom utils packageVersion
#' @import data.table
NULL

.onLoad <- function(libname, pkgname) {
  # Runs when loaded but not attached to search() path; e.g., when a package
  # just Imports (not Depends on) tmle3.
  opts <- list(
    "tmle3.verbose" = FALSE,
    "tmle3.file.path" = tempdir(),
    "tmle3.temp.dir" = tempdir(),
    "tmle3.file.name" = paste0("tmle3-report-", Sys.Date()),
    "tmle3.memoise.learner" = FALSE,
    "tmle3.save.training" = TRUE,
    "tmle3.pcontinuous" = 0.05
  )
  toset <- !(names(opts) %in% names(options()))
  if (any(toset)) options(opts[toset])
  invisible()
}

################################################################################

.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  if (interactive()) {
    v <- utils::packageVersion("tmle3")
    d <- read.dcf(
      system.file("DESCRIPTION", package = "tmle3"),
      fields = c("Packaged", "Built")
    )
    if (is.na(d[1])) {
      if (is.na(d[2])) {
        return() # neither field exists
      } else {
        d <- unlist(strsplit(d[2], split = "; "))[3]
      }
    } else {
      d <- d[1]
    }
    message <- paste0(
      "Please note the package is in early stages of ",
      "development.", "\nCheck often for updates and report",
      "bugs at http://github.com/tlverse/tmle3.", "\n"
    )
    packageStartupMessage("tmle3 ", v)
    packageStartupMessage(message)
  }
}
