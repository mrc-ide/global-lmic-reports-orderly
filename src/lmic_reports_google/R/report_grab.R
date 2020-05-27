reports_day <- function(date = NULL) {
  
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
    message(sprintf("No 'ecdc' report for '%s'", as.character(date)))
    return(NULL)
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
             AND report = "lmic_reports_google"
             AND parameters.name = "iso3c"
           ORDER BY country, report_version.id'
  sql <- sprintf(sql, paste(sprintf('"%s"', id), collapse = ", "))
  reports <- DBI::dbGetQuery(db, sql)
  
  if(nrow(reports) == 0) {
    DBI::dbDisconnect(db)
    return(NULL)
  } else {
  
  if (any(duplicated(reports$country))) {
    keep <- tapply(seq_len(nrow(reports)), reports$country, max)
    reports <- reports[keep, ]
    rownames(reports) <- NULL
  }
  
  reports$date <- as.character(date)
  DBI::dbDisconnect(db)
  return(reports)
  
  }
}


generate_draws <- function(scan_results, squire_model, replicates, n_particles, forecast,
                           country, population, interventions, data) {
  
  # carry out sims drawn from the grid
  if (is.null(scan_results$z)) {
    res <- squire:::sample_grid_scan(scan_results = scan_results,
                                     n_sample_pairs = replicates,
                                     n_particles = n_particles,
                                     forecast_days = forecast ,
                                     full_output = TRUE)
  } else {
    res <- squire:::sample_3d_grid_scan(scan_results = scan_results,
                                        n_sample_pairs = replicates,
                                        n_particles = n_particles,
                                        forecast_days = forecast ,
                                        full_output = TRUE)
  }
  
  
  # recreate model output for each type of model(ish)
  if (inherits(squire_model, "stochastic")) {
    
    # create a fake run object and fill in the required elements
    r <- squire_model$run_func(country = country,
                               contact_matrix_set = scan_results$inputs$model_params$contact_matrix_set,
                               tt_contact_matrix = scan_results$inputs$model_params$tt_matrix,
                               hosp_bed_capacity = scan_results$inputs$model_params$hosp_bed_capacity,
                               tt_hosp_beds = scan_results$inputs$model_params$tt_hosp_beds,
                               ICU_bed_capacity = scan_results$inputs$model_params$ICU_bed_capacity,
                               tt_ICU_beds = scan_results$inputs$model_params$tt_ICU_beds,
                               population = population,
                               replicates = 1,
                               time_period = nrow(res$trajectories))
    
    # first let's create the output
    names(res)[names(res) == "trajectories"] <- "output"
    dimnames(res$output) <- list(dimnames(res$output)[[1]], dimnames(r$output)[[2]], NULL)
    r$output <- res$output
    
    # and adjust the time as before
    full_row <- match(0, apply(r$output[,"time",],2,function(x) { sum(is.na(x)) }))
    saved_full <- r$output[,"time",full_row]
    for(i in seq_len(replicates)) {
      na_pos <- which(is.na(r$output[,"time",i]))
      full_to_place <- saved_full - which(rownames(r$output) == as.Date(max(data$date))) + 1L
      if(length(na_pos) > 0) {
        full_to_place[na_pos] <- NA
      }
      r$output[,"time",i] <- full_to_place
    }
    
  } else if (inherits(squire_model, "deterministic")) {
    r <- list("output" = res$trajectories)
    r <- structure(r, class = "squire_simulation")
  }
  
  # second let's recreate the output
  r$model <- res$inputs$model$odin_model(
    user = res$inputs$model_params, unused_user_action = "ignore"
  )
  
  # we will add the interventions here so that we now what times are needed for projection
  r$interventions <- interventions
  
  # as well as adding the scan_results so it's easy to draw from the scan again in the future
  r$scan_results <- scan_results
  
  # and add the parameters that changed between each simulation, i.e. drawn from gris
  r$replicate_parameters <- res$param_grid
  
  # and fix the replicates
  r$parameters$replicates <- replicates
  r$parameters$time_period <- as.numeric(diff(as.Date(range(rownames(r$output)))))
  
  return(r)
  
}

named_list <- function(...) {
  get <- as.character(match.call())
  l <- list(...)
  names(l) <- tail(get, -1)
  return(l)
}
