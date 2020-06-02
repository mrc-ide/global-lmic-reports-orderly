#!/usr/bin/env Rscript

file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}

## Possibly useful:
## dir("gh-pages", pattern = "index\\.html$", recursive = TRUE)
copy_results <- function(date = NULL, is_latest = TRUE) {
  
  # first set up results dir
  target <- "gh-results/analysis/data/raw_data/server_results"
  #unlink("gh-results", recursive = TRUE, force = TRUE)
  
  # copy over orderly utilities
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  file_copy("orderly.sqlite", target)
  file_copy("orderly_config.yml", target)
  file_copy("global/", target)
  
  # now copy over ecdc, brt and reports
  system("echo pre-DB")
  db <- orderly::orderly_db("destination")
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }
  
  system("echo pre-ecdc")
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
  id_ecdc_max <- max(id)
  
  # copy lmic_reports_google
  src <- file.path("archive", "ecdc", id_ecdc_max)
  dest <- file.path(target,sprintf("%s/%s/%s", "archive", "ecdc",id_ecdc_max))
  worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
  worked <- mapply(file_copy, from = src, to = dirname(dest))
  
  
  ##  --------------------------------------------------------------------------
  ## LMIC REPORTS COPY ---------------------------------------------------------
  ##  --------------------------------------------------------------------------
  
  system("echo pre-report")
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
             AND report = "lmic_reports_google"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql)
  
  if (any(duplicated(reports$country))) {
    keep <- tapply(seq_len(nrow(reports)), reports$country, max)
    reports <- reports[keep, ]
    rownames(reports) <- NULL
  }
  
  reports$date <- as.character(date)
  
  # copy lmic_reports_google key bits
  src <- file.path("archive", "lmic_reports_google", reports$id)
  dest <- file.path(target,sprintf("%s/%s/%s", "archive", "lmic_reports_google",reports$id))
  worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
  
  worked <- mapply(function(from, to){
    append <- c("fitting.pdf", "grid_out.rds", "projections.csv", "orderly_run.rds")
    file_copy(file.path(from, append), to)
  }, from = src, to = dest)
 
  ##  --------------------------------------------------------------------------
  ## BRT COPY ------------------------------------------------------------------
  ##  --------------------------------------------------------------------------
 
  system("echo pre-brt")
  sql <- 'SELECT report_version.id
            FROM report_version
            JOIN parameters
              ON parameters.report_version = report_version.id
           WHERE report_version.report = "brt_google_mobility"
             AND parameters.value = $1'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql, date)
  brt_id_max <- max(reports$id)
  
  # copy lmic_reports_google
  src <- file.path("archive", "brt_google_mobility", brt_id_max)
  dest <- file.path(target,sprintf("%s/%s/%s", "archive", "brt_google_mobility",reports$id))
  worked <- vapply(dest, dir.create, logical(1), recursive = TRUE)
  worked <- mapply(file_copy, from = src, to = dirname(dest))

}


if (!interactive()) {
  usage <- "Usage:\n./copy_results.R [<date>]"
  args <- docopt::docopt(usage)
  copy_results(args$date)
}
