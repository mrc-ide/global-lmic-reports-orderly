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
  ecdc <- readRDS(file.path("archive/ecdc/", max(id), "ecdc_all.rds"))
  
  
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
             AND report = "lmic_reports_google_pmcmc_no_decouple"
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
  pars_init <- readRDS("src/lmic_reports_google_pmcmc_no_decouple/pars_init.rds")
  
  initial_conditions <- for(x in seq_along(reports$id)) {
    
    out <- readRDS(file.path("archive/lmic_reports_google_pmcmc_no_decouple",reports$id[x],"grid_out.rds"))
    mc <- do.call(rbind, lapply(out$pmcmc_results$chains, "[[", "results"))
    best <- mc[which.max(mc$log_posterior),]  
    best$start_date <- as.character(squire:::offset_to_start_date(out$pmcmc_results$inputs$data$date[1], round(best$start_date)))
    best <- best[,1:4]
    best$pld <- out$interventions$date_Meff_change
    rownames(best) <- NULL
    best$iso3c <- reports$country[x]
    
    if(reports$country[x] %in% pars_init$iso3c) {
      pars_init[which(pars_init$iso3c == reports$country[x]),] <- best
    } 
    
  }
  
  saveRDS(pars_init, "src/lmic_reports_google_pmcmc_no_decouple/pars_init.rds")
  
  ## Remove HICs
  rl <- readLines(file.path(here::here(),"countries"))
  to_remove <- stringr::str_sub(rl[(grep("Other HICs", rl) + 1) : length(rl)], -3)
  reports <- reports[!reports$country %in% to_remove,]
  
  ## ---------------------------------------------------------------------------
  ## copy reports --------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  target <- "gh-pages"
  
  src <- file.path("archive", "lmic_reports_google_pmcmc_no_decouple", reports$id)
  dest <- sprintf("gh-pages/%s/%s", reports$country, reports$date)
  copy <- c("index.html",
            "index_files/figure-html",
            "projections.csv",
            "index.pdf",
            "input_params.json")
  
  for (i in seq_along(dest)) {
    message(sprintf("Copying %s (%s)", dest[[i]], reports$id[[i]]))
    dir.create(dest[[i]], FALSE, TRUE)
    file_copy(file.path(src[[i]], copy), dest[[i]])
    if (is_latest) {
      dest_latest <- dirname(dest[[i]])
      prev <- dir(dest_latest, full.names = TRUE, pattern = "\\.")
      unlink(c(prev[-grep("v1|v2", prev)], file.path(dest_latest, "figures")), recursive = TRUE)
      file_copy(dir(dest[[i]], full.names = TRUE), dest_latest)
      
      # remove report if no deaths in last 20 days
      if(sum(head(ecdc[which(ecdc$countryterritoryCode == reports$country[i]),]$deaths,20), na.rm = TRUE)==0) {
        prev <- dir(dest_latest, full.names = TRUE, pattern = "\\.")
        unlink(grep("index", prev, value = TRUE), recursive = TRUE)
      }
      
    }
  }
  
  pdf_input <- file.path(src, "index.pdf")
  message(sprintf("Building combined pdf from %d files", length(pdf_input)))
  qpdf::pdf_combine(pdf_input, "gh-pages/combined_reports.pdf")
  
  ## Aha, this is so naughty, but probably a reasonable shout given
  ## the situation.  The alternative is to depend on _all_ the country
  ## tasks for that date.
  summaries <- do.call(rbind,
                       lapply(file.path(src, "summary_df.rds"), readRDS))
  saveRDS(summaries, "src/index_page/summaries.rds")
  saveRDS(summaries, "src/regional_page/summaries.rds")
  
  projections <- do.call(rbind,
                         lapply(file.path(src, "projections.csv"), read.csv))
  dir.create("gh-pages/data", FALSE, TRUE)
  write.csv(projections, paste0("gh-pages/data/",date,"_v3.csv"), row.names = FALSE, quote = FALSE)
  cwd <- getwd()
  setwd("gh-pages/data/")
  zip(paste0(date,"_v3.csv.zip"),paste0(date,"_v3.csv"))
  file.remove(paste0(date,"_v3.csv"))
  setwd(cwd)
  saveRDS(projections, paste0("src/index_page/all_data.rds"))
  saveRDS(projections, paste0("src/regional_page/all_data.rds"))
  
  ## ---------------------------------------------------------------------------
  ## rt grab -------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  rt <- lapply(seq_along(reports$id), function(x) {
    
    iso <- reports$country[x]
    out <- file.path("archive", "lmic_reports_google_pmcmc_no_decouple", reports$id[x], "grid_out.rds")
    out <- readRDS(out)
    
    # create the Rt data frame
    rts <- lapply(seq_len(length(out$replicate_parameters$R0)), function(y) {
      
      tt <- squire:::intervention_dates_for_odin(dates = out$interventions$date_R0_change, 
                                                 change = out$interventions$R0_change, 
                                                 start_date = out$replicate_parameters$start_date[y],
                                                 steps_per_day = 1/out$parameters$dt)
      
      Rt <- squire:::evaluate_Rt_pmcmc(
        R0_change = tt$change, 
        date_R0_change = tt$dates, 
        R0 = out$replicate_parameters$R0[y], 
        pars = list(
          Meff = out$replicate_parameters$Meff[y],
          Meff_pl = out$replicate_parameters$Meff_pl[y]
        ),
        Rt_args = out$pmcmc_results$inputs$Rt_args) 
      
      df <- data.frame(
        "Rt" = Rt,
        "date" = tt$dates,
        "iso" = iso,
        rep = y,
        stringsAsFactors = FALSE)
      df$pos <- seq_len(nrow(df))
      return(df)
    } )
    
    rt <- do.call(rbind, rts)
    return(rt)
  })
  names(rt) <- reports$country
  
  rt_all <- do.call(rbind, rt)
  rt_all$date <- as.Date(rt_all$date)
  rt_all <- rt_all[,c(3,2,1,4,5)]
  
  library(magrittr)
  date_0 <- as.Date(date)
  new_rt_all <- rt_all %>%
    dplyr::group_by(iso, rep) %>% 
    dplyr::arrange(date) %>% 
    tidyr::complete(date = seq.Date(min(rt_all$date), date_0, by = "days")) 
  
  column_names <- colnames(new_rt_all)[-c(1,2,3)]
  new_rt_all <- tidyr::fill(new_rt_all, tidyselect::all_of(column_names), .direction = c("down"))
  new_rt_all <- tidyr::fill(new_rt_all, tidyselect::all_of(column_names), .direction = c("up"))
  
  sum_rt <- dplyr::group_by(new_rt_all, iso, date) %>% 
    dplyr::summarise(Rt_min = quantile(Rt, 0.025),
                     Rt_q25 = quantile(Rt, 0.25),
                     Rt_q75 = quantile(Rt, 0.75),
                     Rt_max = quantile(Rt, 0.975),
                     Rt = median(Rt)) 
  sum_rt$continent <- countrycode::countrycode(sum_rt$iso, "iso3c", "continent")
  saveRDS(sum_rt, paste0("src/index_page/sum_rt.rds"))
  saveRDS(sum_rt, paste0("src/regional_page/sum_rt.rds"))
  
  
  
}


if (!interactive()) {
  usage <- "Usage:\n./copy_outputs.R [<date>]"
  args <- docopt::docopt(usage)
  copy_outputs(args$date)
}
