update_gh_fits <- function(ids, filename){
  repo <- here::here()
  destination <- file.path(
    repo, "gh-fits",
    "standard"
  )

  #upload to folder replacing existing files
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)

  #minimise RAM useage by loading and saving in turn
  purrr::walk(ids, function(x){
    fit <- readRDS(file.path(repo, "archive", task, x$id, filename))
    saveRDS(fit,
            file.path(
              destination,
              paste0(x$iso3c, ".Rds")
            ))
  })
}
