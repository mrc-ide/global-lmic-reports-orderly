update_gh_fits <- function(ids, task) {
  repo <- here::here()
  destination_excess <- file.path(
    repo, "gh-fits"
  )
  excess_name <- "excess_out.Rds"
  #upload to folder replacing existing files
  dir.create(destination_excess, recursive = TRUE, showWarnings = FALSE)

  #minimise RAM useage by loading and saving in turn
  purrr::walk(ids, function(x) {
    tryCatch(readRDS(file.path(repo, "archive", task, x$id, excess_name)) %>%
    saveRDS(file.path(
      destination_excess,
              paste0(x$iso3c, ".Rds")
    )), error = function(e) NULL)
    }
  )
}
