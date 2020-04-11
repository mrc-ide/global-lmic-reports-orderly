#!/usr/bin/env Rscript

## TODO:
## * overall index page
## * ISO3 country codes
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
copy_outputs <- function(date = Sys.Date(), is_latest = TRUE) {
  db <- orderly::orderly_db("destination")

  sql <- 'SELECT report_version.id, report_version.date, parameters.value as country
  FROM report_version
  JOIN parameters
    ON parameters.report_version = report_version.id
 WHERE report_version.report = "lmic_reports"'

  reports <- DBI::dbGetQuery(db, sql)
  reports <- reports[reports$date >= date & reports$date < date + 1, ]

  if (any(duplicated(reports$country))) {
    stop("OJ needs to filter this by most recent per country")
  }

  ## Remove this once the lmic task takes ISO3
  remap <- c(Angola = "AGO", Senegal = "SEN")
  reports$country <- unname(remap[reports$country])

  target <- "gh-pages"

  src <- file.path("archive", "lmic_reports", reports$id)
  dest <- sprintf("gh-pages/%s/%s", reports$country, date)
  copy <- c("report.html",
            # "figures" # uncomment once it exists
            # "fig1.png", "fig2.png",
            "report.pdf")

  for (i in seq_along(dest)) {
    dir.create(dest, FALSE, TRUE)
    file_copy(file.path(src[[i]], copy), dest[[i]])
    ## Remove after report.html renamed to index.html
    file.rename(file.path(dest[[i]], "report.html"),
                file.path(dest[[i]], "index.html"))
    if (is_latest) {
      dest_latest <- dirname(dest)
      prev <- dir(dest_latest, pattern = "\\.")
      unlink(c(prev, file.path(dest_latest, "figures")), recursive = TRUE)
      file_copy(dir(dest[[i]], full.names = TRUE), dest_latest)
    }
  }
}


if (!interactive()) {
  copy_outputs()
}
