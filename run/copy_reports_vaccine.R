#!/usr/bin/env Rscript

file_copy <- function(from, to) {
  ok <- file.copy(from, to, overwrite = TRUE, recursive = TRUE)
  if (any(!ok)) {
    stop("There was an error copying files")
  }
}

get_immunity_ratios <- function(out) {
  
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$pmcmc_results$inputs$model_params$contact_matrix_set[[1]],
    out$pmcmc_results$inputs$model_params$population
  )
  
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  
  # assertions
  squire:::assert_single_pos(dur_ICase, zero_allowed = FALSE)
  squire:::assert_single_pos(dur_IMild, zero_allowed = FALSE)
  squire:::assert_numeric(prob_hosp)
  squire:::assert_numeric(mixing_matrix)
  squire:::assert_square_matrix(mixing_matrix)
  squire:::assert_same_length(mixing_matrix[,1], prob_hosp)
  
  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }
  
  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }
  
  index <- squire:::odin_index(out$model)
  pop <- out$parameters$population
  t_now <- which(as.Date(rownames(out$output)) == max(out$pmcmc_results$inputs$data$date))
  prop_susc <- lapply(seq_len(dim(out$output)[3]), function(x) {
    t(t(out$output[seq_len(t_now), index$S, x])/pop)
  } )
  
  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  
  adjusted_eigens <- lapply(prop_susc, function(x) {
    
    unlist(lapply(seq_len(nrow(x)), function(y) {
      if(any(is.na(x[y,]))) {
        return(NA)
      } else {
        Re(eigen(mixing_matrix*x[y,]*relative_R0_by_age)$values[1])
      }
    }))
    
  })
  
  betas <- lapply(out$replicate_parameters$R0, function(x) {
    squire:::beta_est(squire_model = out$pmcmc_results$inputs$squire_model, 
                      model_params = out$pmcmc_results$inputs$model_params, 
                      R0 = x)
  })
  
  ratios <- lapply(seq_along(betas), function(x) {
    (betas[[x]] * adjusted_eigens[[x]]) / out$replicate_parameters$R0[[x]]
  })
  
  return(ratios)
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
             AND report = "lmic_reports_vaccine"
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
  pars_init <- readRDS("src/lmic_reports_vaccine/pars_init.rds")
  pars <- vector("list", length(reports$id))
  also_to_go <- vector("list", length(reports$id))
  for(x in seq_along(reports$id)) {
    
    out <- readRDS(file.path("archive/lmic_reports_vaccine",reports$id[x],"grid_out.rds"))
    if("pmcmc_results" %in% names(out)) {
    mc <- do.call(rbind, lapply(out$pmcmc_results$chains, "[[", "results"))
    best <- mc[which.max(mc$log_posterior),]  
    best <- best[,seq_len(ncol(best)-3)]
    rownames(best) <- NULL
    best$start_date <- as.character(squire:::offset_to_start_date(out$pmcmc_results$inputs$data$date[1], round(best$start_date)))
    best$iso3c <- reports$country[x]
    best$date_Meff_change <- out$pmcmc_results$inputs$Rt_args$date_Meff_change
    best$Rt_shift_duration <- out$pmcmc_results$inputs$Rt_args$Rt_shift_duration
    best$Rt_rw_duration <- out$pmcmc_results$inputs$Rt_args$Rt_rw_duration
    best$covariance_matrix <- out$pmcmc_results$chains$chain1$covariance_matrix[1]
    best$scaling_factor <- mean(
      c(tail(na.omit(out$pmcmc_results$chains$chain1$scaling_factor),1),
        tail(na.omit(out$pmcmc_results$chains$chain2$scaling_factor),1),
        tail(na.omit(out$pmcmc_results$chains$chain3$scaling_factor),1))
      )
    
    # for now combine here
    pars[[x]] <- best
    also_to_go[[x]] <- NULL
    } else {
      pars[[x]] <- NULL
      also_to_go[[x]] <- reports$country[[x]]
    }

  }
  
  #names(pars)[unlist(lapply(lapply(pars,is.null), isFALSE))] <- reports$country[[unlist(lapply(lapply(pars,is.null), isFALSE))]]
  names(pars) <- reports$country
  
  saveRDS(pars, "src/lmic_reports_vaccine/pars_init.rds")
  
  ## Remove HICs
  rl <- readLines(file.path(here::here(),"countries"))
  to_remove <- stringr::str_sub(rl[(grep("Other HICs", rl) + 1) : length(rl)], -3)
  to_remove <- c(to_remove, unlist(also_to_go))
  hic_pos <- which(reports$country %in% to_remove)
  
  ## ---------------------------------------------------------------------------
  ## copy reports --------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  target <- "gh-pages"
  
  src <- file.path("archive", "lmic_reports_vaccine", reports$id)
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
  message(sprintf("Building combined input params from %d files", length(all_jsons)))
  zip("gh-pages/combined_input_params_jsons.zip", files = all_jsons)
  
  # from here remove hic
  pdf_input <- file.path(src[-hic_pos], "index.pdf")
  message(sprintf("Building combined pdf from %d files", length(pdf_input)))
  qpdf::pdf_combine(pdf_input, "gh-pages/combined_reports.pdf")
  
  ## Aha, this is so naughty, but probably a reasonable shout given
  ## the situation.  The alternative is to depend on _all_ the country
  ## tasks for that date.
  summaries <- do.call(rbind,
                       lapply(file.path(src[-hic_pos], "summary_df.rds"), readRDS))
  saveRDS(summaries, "src/index_page/summaries.rds")
  saveRDS(summaries, "src/regional_page/summaries.rds")
  
  projections <- do.call(rbind,
                         lapply(file.path(src, "projections.csv"), read.csv))
  dir.create("gh-pages/data", FALSE, TRUE)
  projections$version <- "v8"
  write.csv(projections, paste0("gh-pages/data/",date,"_v8.csv"), row.names = FALSE, quote = FALSE)
  cwd <- getwd()
  setwd("gh-pages/data/")
  zip(paste0(date,"_v8.csv.zip"),paste0(date,"_v8.csv"))
  file.remove(paste0(date,"_v8.csv"))
  setwd(cwd)
  
  hic_pos_projections <- which(projections$iso3c %in% to_remove)
  
  saveRDS(projections[-hic_pos_projections,], paste0("src/index_page/all_data.rds"))
  saveRDS(projections[-hic_pos_projections,], paste0("src/regional_page/all_data.rds"))
  
  
}


if (!interactive()) {
  usage <- "Usage:\n./copy_outputs.R [<date>]"
  args <- docopt::docopt(usage)
  copy_outputs(args$date)
}
