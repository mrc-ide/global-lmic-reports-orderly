update_gh_fits <- function(ids){
  repo <- here::here()
  destination_excess <- file.path(
    repo, "gh-fits",
    "excess_mortality"
  )
  excess_name <- "excess_deaths.Rds"
  destination_reported <- file.path(
    repo, "gh-fits",
    "reported_deaths"
  )
  reported_name <- "excess_deaths.Rds"
  #upload to folder replacing existing files
  dir.create(c(destination_excess, destination_reported), recursive = TRUE, showWarnings = FALSE)

  #minimise RAM useage by loading and saving in turn
  purrr::walk(ids, function(x){
    silent <- readRDS(file.path(repo, "archive", task, x$id, excess_name)) %>%
    saveRDS(file.path(
      destination_excess,
              paste0(x$iso3c, ".Rds")
            ))
    readRDS(file.path(repo, "archive", task, x$id, reported_name)) %>%
      saveRDS(file.path(
        destination_reported,
        paste0(x$iso3c, ".Rds")
      ))
  })
}
