update_gh_fits <- function(ids){
  repo <- here::here()
  destination_excess <- file.path(
    repo, "gh-fits",
    "excess_mortality"
  )
  excess_name <- "excess_out.Rds"
  destination_reported <- file.path(
    repo, "gh-fits",
    "reported_deaths"
  )
  reported_name <- "reported_out.Rds"
  #upload to folder replacing existing files
  dir.create(destination_excess, recursive = TRUE, showWarnings = FALSE)
  dir.create(destination_reported, recursive = TRUE, showWarnings = FALSE)

  #minimise RAM useage by loading and saving in turn
  purrr::walk(ids, function(x){
    silent <- tryCatch(readRDS(file.path(repo, "archive", task, x$id, excess_name)) %>%
    saveRDS(file.path(
      destination_excess,
              paste0(x$iso3c, ".Rds")
    )), error = function(e){NULL}
    )
  tryCatch(
    readRDS(file.path(repo, "archive", task, x$id, reported_name)) %>%
      saveRDS(file.path(
        destination_reported,
        paste0(x$iso3c, ".Rds")
      )), error = function(e){NULL}
  )
  })
}
