#!/usr/bin/env Rscript

fits_for_checking <- function(date, report = "lmic_reports_vaccine") {

# ------------------------------------------------------------------------------
### Get the reports ------------------------------------------------------------
# ------------------------------------------------------------------------------
reports_day <- function(date = NULL, report = "lmic_reports_vaccine") {

  db <- orderly::orderly_db("destination")
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }

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
  sql <- paste0('SELECT report_version.id, parameters.value as country
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
             AND report = "', report, '"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id')
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql)

  if (any(duplicated(reports$country))) {
    keep <- tapply(seq_len(nrow(reports)), reports$country, max)
    reports <- reports[keep, ]
    rownames(reports) <- NULL
  }

  reports$date <- as.character(date)

  return(reports)
}

reports <- reports_day(date, report)
if(grepl("vaccine", report)){
  file_suffix <- "_vaccine"
} else if(is.na(strsplit(report,"_")[[1]][3])){
  file_suffix <- ""
} else if(is.na(strsplit(report,"_")[[1]][4])){
  file_suffix <- "_google"
} else if(is.na(strsplit(report,"_")[[1]][5])){
  file_suffix <- "_pmcmc"
} else if(is.na(strsplit(report,"_")[[1]][6])){
  file_suffix <- "_pmcmc_6p"
}  else {
  file_suffix <- "_spline"
}

pdfs <- file.path(paste0("archive/",report), reports$id, "fitting.pdf")
fz <- file.size(pdfs)
dir.create("fits", showWarnings = FALSE)
td <- tempdir()
dir.create(file.path(td, "fits"))
pdf_out <- file.path(td, "fits",paste0(reports$country[which(fz>0)],file_suffix,".pdf"))
pdf_input <- file.copy(pdfs[which(fz>0)], pdf_out, overwrite = TRUE)
qpdf::pdf_combine(pdf_out, file.path(paste0("fits/",report,"_",date,".pdf")))

}

if (!interactive()) {
  usage <- "Usage:\n./fits_for_checking.R [<date>] [<report>]"
  args <- docopt::docopt(usage)
  fits_for_checking(args$date, args$report)
}
