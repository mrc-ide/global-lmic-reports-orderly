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
  date <- as.character(date)

  ## TODO: - use date only

  ## First find the id corresponding to the ecdc report with data.  If
  ## there are more than one, it's not totally clear what you want to
  ## do as you might want to take the earliest or the latest.
  ## Probably we want to take *all* and do the join over that, which
  ## is easy enough to do if you replace the '= $1' and replace with
  ## 'IN (%s)' and interpolate 'paste(sprintf('"%s"', id), collapse = ", ")'
  sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "ecdc"
             AND parameters.value = $1'
  id <- DBI::dbGetQuery(db, sql, date)$id
  if (length(id) == 0L) {
    stop(sprintf("No 'ecdc' report for '%s'", as.character(date)))
  } else if (length(id) > 1) {
    message(sprintf("Multiple 'ecdc' reports for '%s'", as.character(date)))
  }

  ## Then find all lmic_reports reports that use files from this ecdc
  ## report.  This is a bit awful and I might add direct link or a
  ## view to make this easier at some point.
  sql <- 'SELECT report_version.id, parameters.value as country
            FROM report_version_artefact
            JOIN file_artefact
              ON file_artefact.artefact = report_version_artefact.id
            JOIN depends
              ON depends.use = file_artefact.id
            JOIN report_version
              ON report_version.id = depends.report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version_artefact.report_version IN (%s)
             AND report = "lmic_reports"
           ORDER BY country, report_version.id'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql)

  if (any(duplicated(reports$country))) {
    keep <- tapply(seq_len(nrow(reports)), reports$country, max)
    reports <- reports[keep, ]
    rownames(reports) <- NULL
  }

  ## Remove this once the lmic task takes ISO3
  remap <- c(Angola = "AGO", Senegal = "SEN")
  reports$country <- unname(remap[reports$country])
  reports$date <- as.character(date)

  target <- "gh-pages"

  src <- file.path("archive", "lmic_reports", reports$id)
  dest <- sprintf("gh-pages/%s/%s", reports$country, reports$date)
  copy <- c("report.html",
            # "figures" # uncomment once it exists
            # "fig1.png", "fig2.png",
            "report.pdf")

  for (i in seq_along(dest)) {
    dir.create(dest[[i]], FALSE, TRUE)
    file_copy(file.path(src[[i]], copy), dest[[i]])
    ## Remove after report.html renamed to index.html
    file.rename(file.path(dest[[i]], "report.html"),
                file.path(dest[[i]], "index.html"))
    if (is_latest) {
      dest_latest <- dirname(dest[[i]])
      prev <- dir(dest_latest, pattern = "\\.")
      unlink(c(prev, file.path(dest_latest, "figures")), recursive = TRUE)
      file_copy(dir(dest[[i]], full.names = TRUE), dest_latest)
    }
  }
}


if (!interactive()) {
  copy_outputs()
}
