#!/usr/bin/env Rscript

## TODO:
## * overall index page
## * figures updated and change below
## * filter by most recent run per country

file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}

## Possibly useful:
## dir("gh-pages", pattern = "index\\.html$", recursive = TRUE)
copy_index <- function(date = Sys.Date(), is_latest = TRUE) {
  
  db <- orderly::orderly_db("destination")
  index <- orderly:::get_ids_by_name(db, "index_page")
  params <- orderly:::get_ids_by_name(db, "parameters")
  
  target <- "gh-pages"
  
  src_index <- file.path("archive", "index_page", index$id[which.max(as.Date(index$date))],"index.html")
  src_params <- file.path("archive", "parameters", params$id[which.max(as.Date(params$date))],"parameters.html")
  src <- c(src_index, src_params)
  dest <- c(sprintf("gh-pages/%s","index.html"),
            sprintf("gh-pages/%s","parameters.html") )
            
    file_copy(file.path(src), dest)
  
}


if (!interactive()) {
  copy_index()
}

