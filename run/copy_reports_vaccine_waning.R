#!/usr/bin/env Rscript

file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}


## Possibly useful:
## dir("gh-pages", pattern = "index\\.html$", recursive = TRUE)
copy_outputs <- function(date = NULL, is_latest = TRUE) {
  db <- orderly::orderly_db("destination")
  if (is.null(date)) {
    date <- as.character(Sys.Date())
  }

  ## ---------------------------------------------------------------------------
  ## reports grab --------------------------------------------------------------
  ## ---------------------------------------------------------------------------

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
  ecdc <- readRDS(file.path("archive/ecdc/", max(id), "jhu_all.rds"))


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
             AND report = "lmic_reports_vaccine_booster"
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

  ## ---------------------------------------------------------------------------
  ## initial conditions --------------------------------------------------------
  ## ---------------------------------------------------------------------------

  # get old conditions
  also_to_go <- vector("list", length(reports$id))
  for(x in seq_along(reports$id)) {
    out <- readRDS(file.path("archive/lmic_reports_vaccine_booster",reports$id[x],"grid_out.rds"))
    if("pmcmc_results" %in% names(out)) {
      also_to_go[[x]] <- NULL
    } else {
      also_to_go[[x]] <- reports$country[[x]]
    }

  }


  ## Remove HICs
  rl <- readLines(file.path(here::here(),"countries"))
  to_remove <- stringr::str_sub(rl[(grep("Other HICs", rl) + 1) : length(rl)], -3)
  to_remove <- c(to_remove, unlist(also_to_go))
  non_hic_pos <- which(!(reports$country %in% to_remove))

  ## ---------------------------------------------------------------------------
  ## copy reports --------------------------------------------------------------
  ## ---------------------------------------------------------------------------

  target <- "gh-pages"

  src <- file.path("archive", "lmic_reports_vaccine_booster", reports$id)
  dest <- sprintf("gh-pages/%s/%s", reports$country, reports$date)
  copy <- c("index.html",
            "projections.csv",
            "index.pdf",
            "input_params.json")

  for (i in seq_along(dest)) {
    message(sprintf("Copying %s (%s)", dest[[i]], reports$id[[i]]))
    dir.create(dest[[i]], FALSE, TRUE)
    to_copy <- file.path(src[[i]], copy)
    fz <- file.size(to_copy)
    to_copy <- to_copy[fz>0]
    file_copy(to_copy, dest[[i]])
    if (is_latest) {
      dest_latest <- dirname(dest[[i]])
      prev <- dir(dest_latest, full.names = TRUE, pattern = "\\.")
      unlink(c(prev, file.path(dest_latest, "figures")), recursive = TRUE)
      file_copy(dir(dest[[i]], full.names = TRUE), dest_latest)

      # remove report if no deaths in last 20 days
      if(sum(tail(ecdc[which(ecdc$countryterritoryCode == reports$country[i]),]$deaths,20), na.rm = TRUE)==0) {
        prev <- dir(dest_latest, full.names = TRUE, pattern = "\\.")
        unlink(grep("index", prev, value = TRUE), recursive = TRUE)
      }

    }
  }

  # zip all inputs together for Rich
  all_jsons <- file.path(src, "input_params.json")
  json_tos <- file.path(paste0(reports$country, "_", "input_params.json"))
  file.copy(all_jsons, json_tos)
  message(sprintf("Building combined input params from %d files", length(all_jsons)))
  dir.create("gh-pages/input_jsons", showWarnings = FALSE)
  zip(zipfile = paste0("gh-pages/input_jsons/",date, "_combined_input_params_jsons.zip"), files = json_tos)
  file.remove(json_tos)

  # from here remove hic as the website is only really meant to be be LMIC countries
  pdf_input <- file.path(src[non_hic_pos], "index.pdf")
  message(sprintf("Building combined pdf from %d files", length(pdf_input)))
  qpdf::pdf_combine(pdf_input, "gh-pages/combined_reports.pdf")

  ## Aha, this is so naughty, but probably a reasonable shout given
  ## the situation.  The alternative is to depend on _all_ the country
  ## tasks for that date.
  summaries <- do.call(rbind,
                       lapply(file.path(src[non_hic_pos], "summary_df.rds"), readRDS))
  saveRDS(summaries, "src/index_page/summaries.rds")
  saveRDS(summaries, "src/regional_page/summaries.rds")

  projections <- do.call(rbind,
                         lapply(file.path(src, "projections.csv"), read.csv))
  #split into two csv files
  #file with old scenarios
  library(magrittr)
  projections_1 <-
    dplyr::filter(projections,
                  scenario %in% c(
                    "Maintain Status Quo",
                    "Relax Interventions 50%",
                    "Additional 50% Reduction",
                    "Surged Maintain Status Quo",
                    "Surged Additional 50% Reduction",
                    "Surged Relax Interventions 50%"
                  )) %>%
    dplyr::filter(as.Date(date) > min(as.Date(date) + 150)) #TEMP RESOLVE THIS ISSUE
  #new scenarios
  projections_2 <-
    dplyr::filter(projections,
                  scenario %in% c(
                    "Maintain Status Quo",
                    "Optimistic",
                    "Pessimistic",
                    "Surged Maintain Status Quo",
                    "Surged Optimistic",
                    "Surged Pessimistic"
                  )) %>%
    dplyr::filter(as.Date(date) > min(as.Date(date) + 150)) #TEMP RESOLVE THIS ISSUE

  dir.create("gh-pages/data", FALSE, TRUE)
  cwd <- getwd()
  setwd("gh-pages/data/")

  projections_1$version <- "v8"
  write.csv(projections_1, paste0(date,"_v8.csv"), row.names = FALSE, quote = FALSE)
  zip(paste0(date,"_v8.csv.zip"),paste0(date,"_v8.csv"))
  file.remove(paste0(date,"_v8.csv"))

  projections_2$version <- "v9"
  write.csv(projections_2, paste0(date,"_v9.csv"), row.names = FALSE, quote = FALSE)
  zip(paste0(date,"_v9.csv.zip"),paste0(date,"_v9.csv"))
  file.remove(paste0(date,"_v9.csv"))

  setwd(cwd)

  all_data <-  dplyr::filter(projections,
                             scenario %in% c(
                               "Maintain Status Quo",
                               "Optimistic",
                               "Pessimistic",
                               "Surged Maintain Status Quo",
                               "Surged Optimistic",
                               "Surged Pessimistic"
                             )) %>%
    dplyr::mutate(version = "v9") %>%
    dplyr::filter(!iso3c %in% to_remove)

  saveRDS(all_data, paste0("src/index_page/all_data.rds"))
  saveRDS(all_data, paste0("src/regional_page/all_data.rds"))


}


if (!interactive()) {
  usage <- "Usage:\n./copy_reports_vaccine_waning.R [<date>]"
  args <- docopt::docopt(usage)
  copy_outputs(args$date)
}
